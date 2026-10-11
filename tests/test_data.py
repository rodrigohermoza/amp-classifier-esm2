"""Tests del parseo y la limpieza.

Varios son regresiones de problemas reales documentados en
docs/09-problemas-encontrados.md.
"""

import pandas as pd

from ampesm.data.clean import clean_records, deduplicate, is_standard, match_lengths, seq_sha1
from ampesm.data.parse import PeptideRecord, iter_fasta, parse_fasta_source
from ampesm.data.sources import get_sources

APD3 = get_sources(["apd3"])[0]


def rec(seq, source="apd3", label=1, rec_id=None):
    return PeptideRecord(id=rec_id or seq[:6], sequence=seq, source=source, label=label)


class TestFastaParsing:
    def test_basic_fasta(self, tmp_path):
        path = tmp_path / "x.fasta"
        path.write_text(">AP00001\nGIGKFL\n>AP00002\nKKRRKK\n")
        assert iter_fasta(path) == [("AP00001", "GIGKFL"), ("AP00002", "KKRRKK")]

    def test_multiline_sequence_is_joined(self):
        """Un FASTA valido puede partir la secuencia en varias lineas."""
        import tempfile
        from pathlib import Path

        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.fasta"
            path.write_text(">a\nGIGK\nFLHS\nAKKF\n")
            assert iter_fasta(path) == [("a", "GIGKFLHSAKKF")]

    def test_apd3_spurious_header_is_dropped(self, tmp_path):
        """Regresion de P-01.

        El FASTA de APD3 empieza con '>Your search led to 3306 peptides', una
        cabecera sin secuencia. Si no se descarta, desplaza todos los registros
        y envenena el conjunto en silencio.
        """
        path = tmp_path / "apd3.fasta"
        path.write_text(
            ">Your search led to 3306 peptides\n"
            ">AP00001\nGLWSKIKEVGKEAAKAAAKAAGKAALGAVSEAV\n"
            ">AP00002\nYVPLPNVPQPGRRPFPTFPGQGPFNPKIKWPQGY\n"
        )
        records, report = parse_fasta_source(path, APD3)

        assert report["dropped_empty"] == 1
        assert len(records) == 2
        assert records[0].id == "AP00001"
        assert records[0].sequence.startswith("GLWSKIKEV")

    def test_uniprot_header_yields_accession(self, tmp_path):
        path = tmp_path / "u.fasta"
        path.write_text(">sp|P12345|TEST_HUMAN Some protein\nGIGKFL\n")
        records, _ = parse_fasta_source(path, APD3)
        assert records[0].id == "P12345"

    def test_blank_lines_are_ignored(self, tmp_path):
        path = tmp_path / "x.fasta"
        path.write_text(">a\nGIGK\n\n\n>b\nKKRR\n\n")
        assert iter_fasta(path) == [("a", "GIGK"), ("b", "KKRR")]


class TestAlphabetFilter:
    def test_standard_sequence_passes(self):
        assert is_standard("GIGKFLHSAKKFGKAFVGEIMNS")

    def test_unknown_residue_rejected(self):
        assert not is_standard("GIGKXLHS")  # X = desconocido

    def test_ambiguous_residues_rejected(self):
        for bad in ("GIGKBLHS", "GIGKZLHS", "GIGKJLHS"):
            assert not is_standard(bad)

    def test_stop_codon_and_gaps_rejected(self):
        assert not is_standard("GIGK*LHS")
        assert not is_standard("GIGK-LHS")

    def test_empty_rejected(self):
        assert not is_standard("")

    def test_lowercase_rejected(self):
        """La normalizacion a mayusculas es responsabilidad de quien llama."""
        assert not is_standard("gigkflhs")


class TestCleanRecords:
    def test_length_bounds(self):
        records = [rec("AAAA"), rec("AAAAA"), rec("A" * 50), rec("A" * 51)]
        kept, stats = clean_records(records, min_length=5, max_length=50)
        assert len(kept) == 2
        assert stats.counts["apd3"]["dropped_too_short"] == 1
        assert stats.counts["apd3"]["dropped_too_long"] == 1

    def test_non_standard_counted(self):
        kept, stats = clean_records([rec("GIGKXLHSAK")], 5, 50)
        assert kept == []
        assert stats.counts["apd3"]["dropped_non_standard"] == 1

    def test_whitespace_is_stripped(self):
        kept, _ = clean_records([rec("GIG KFL HSA")], 5, 50)
        assert len(kept) == 1
        assert kept[0].sequence == "GIGKFLHSA"

    def test_counts_add_up(self):
        """input = kept + todos los descartes. Cuadrar los conteos es la
        comprobacion mas barata de que nada se pierde en silencio (leccion de P-01)."""
        records = [rec("AAAA"), rec("GIGKFLHSAK"), rec("GIGKXLHSAK"), rec("A" * 60)]
        kept, stats = clean_records(records, 5, 50)
        counts = stats.counts["apd3"]
        dropped = sum(v for k, v in counts.items() if k.startswith("dropped"))
        assert counts["input"] == counts["kept"] + dropped
        assert counts["kept"] == len(kept)


class TestDeduplicate:
    def test_exact_duplicates_collapse(self):
        records = [rec("GIGKFLHSAK", "apd3"), rec("GIGKFLHSAK", "dramp")]
        kept, report = deduplicate(records)
        assert len(kept) == 1
        assert report["duplicate_records_removed"] == 1
        assert kept[0].extra["n_sources"] == 2

    def test_source_priority_respected(self):
        """Ante un empate se conserva la fuente mas curada (APD3 sobre DRAMP)."""
        records = [rec("GIGKFLHSAK", "dramp"), rec("GIGKFLHSAK", "apd3")]
        kept, _ = deduplicate(records)
        assert kept[0].source == "apd3"

    def test_label_conflict_resolves_to_positive(self):
        """Si una secuencia aparece como AMP y como no-AMP, gana el AMP.

        La anotacion ausente en UniProt es el error, no la actividad verificada
        en APD3. El conteo de estos casos mide por lo bajo la tasa de falsos
        negativos del conjunto de negativos.
        """
        records = [
            rec("GIGKFLHSAK", "apd3", label=1),
            rec("GIGKFLHSAK", "uniprot_neg", label=0),
        ]
        kept, report = deduplicate(records)
        assert len(kept) == 1
        assert kept[0].label == 1
        assert report["label_conflicts"] == 1

    def test_distinct_sequences_are_kept(self):
        kept, report = deduplicate([rec("GIGKFLHSAK"), rec("KKRRKKRRKK")])
        assert len(kept) == 2
        assert report["duplicate_records_removed"] == 0

    def test_hash_is_sequence_based(self):
        """El hash depende de la secuencia, no del id: es la clave del cache."""
        assert seq_sha1("GIGKFL") == seq_sha1("GIGKFL")
        assert seq_sha1("GIGKFL") != seq_sha1("GIGKFW")


class TestLengthMatching:
    def _frame(self):
        rows = []
        # 10 positivos de longitud 20, 5 de longitud 30.
        for i in range(10):
            rows.append({"id": f"p20_{i}", "length": 20, "label": 1})
        for i in range(5):
            rows.append({"id": f"p30_{i}", "length": 30, "label": 1})
        # Negativos de sobra en 20 y 30, y 50 en una longitud sin positivos.
        for i in range(100):
            rows.append({"id": f"n20_{i}", "length": 20, "label": 0})
        for i in range(100):
            rows.append({"id": f"n30_{i}", "length": 30, "label": 0})
        for i in range(50):
            rows.append({"id": f"n45_{i}", "length": 45, "label": 0})
        return pd.DataFrame(rows)

    def test_ratio_respected_per_length(self):
        out, report = match_lengths(self._frame(), ratio=3.0, seed=42)
        assert len(out[(out.label == 0) & (out.length == 20)]) == 30
        assert len(out[(out.label == 0) & (out.length == 30)]) == 15
        assert report["negatives_after"] == 45

    def test_lengths_without_positives_are_dropped(self):
        """Negativos de una longitud sin positivos no aportan y desbalancean."""
        out, _ = match_lengths(self._frame(), ratio=3.0, seed=42)
        assert len(out[out.length == 45]) == 0

    def test_positives_are_never_dropped(self):
        out, _ = match_lengths(self._frame(), ratio=3.0, seed=42)
        assert (out.label == 1).sum() == 15

    def test_deficit_is_reported(self):
        """Si falta material no se compensa con otras longitudes; se declara."""
        df = pd.DataFrame(
            [{"id": f"p_{i}", "length": 20, "label": 1} for i in range(10)]
            + [{"id": f"n_{i}", "length": 20, "label": 0} for i in range(5)]
        )
        out, report = match_lengths(df, ratio=3.0, seed=42)
        assert report["deficit"] == 25  # queria 30, habia 5
        assert (out.label == 0).sum() == 5

    def test_deterministic(self):
        a, _ = match_lengths(self._frame(), ratio=3.0, seed=1)
        b, _ = match_lengths(self._frame(), ratio=3.0, seed=1)
        assert a.id.tolist() == b.id.tolist()
