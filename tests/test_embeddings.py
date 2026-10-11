"""Tests de la agregacion de embeddings.

El test que importa es test_padding_does_not_change_embedding. El error que
detecta -- promediar sin enmascarar el relleno -- no produce ningun fallo
visible: simplemente mete la longitud del peptido dentro del embedding, y
entonces el clasificador acierta por el motivo equivocado.

Ver docs/06-embeddings-esm2.md, seccion de agregacion.
"""

import importlib.util

import numpy as np
import pytest

from ampesm.embeddings import MODEL_DIMS, cache_path, load_cache, save_cache

torch_missing = importlib.util.find_spec("torch") is None


class TestCache:
    def test_path_encodes_configuration(self):
        """Cambiar capa o agregacion debe dar otro archivo de cache.

        Si no, una corrida reutilizaria en silencio embeddings calculados con
        otra configuracion y la comparacion entre modelos seria falsa.
        """
        a = cache_path("facebook/esm2_t6_8M_UR50D", "last", "mean")
        b = cache_path("facebook/esm2_t6_8M_UR50D", 3, "mean")
        c = cache_path("facebook/esm2_t6_8M_UR50D", "last", "mean_max")
        d = cache_path("facebook/esm2_t12_35M_UR50D", "last", "mean")
        assert len({a, b, c, d}) == 4

    def test_path_has_no_slashes_from_model_name(self):
        path = cache_path("facebook/esm2_t6_8M_UR50D", "last", "mean")
        assert "/" not in path.name

    def test_roundtrip(self, tmp_path):
        path = tmp_path / "cache.npz"
        original = {"abc123": np.arange(5, dtype=np.float32), "def456": np.ones(5, np.float32)}
        save_cache(path, original)
        loaded = load_cache(path)
        assert set(loaded) == set(original)
        for key in original:
            assert np.allclose(loaded[key], original[key])

    def test_missing_cache_is_empty_not_error(self, tmp_path):
        assert load_cache(tmp_path / "no_existe.npz") == {}

    def test_no_part_file_left_behind(self, tmp_path):
        """La escritura es atomica via .part; no debe quedar basura."""
        path = tmp_path / "cache.npz"
        save_cache(path, {"a": np.ones(3, np.float32)})
        assert path.exists()
        assert not list(tmp_path.glob("*.part"))


@pytest.mark.skipif(torch_missing, reason="PyTorch no instalado")
class TestPooling:
    def test_padding_does_not_change_embedding(self):
        """Dos secuencias en un lote con relleno deben dar lo mismo que solas.

        Este es EL test de esta etapa. Un promedio sin mascara divide por la
        longitud del tensor en vez de por el numero de residuos reales, asi que
        el vector de un peptido corto queda escalado por un factor que depende
        de su longitud. El modelo lee la longitud y el resultado parece bueno.
        """
        import torch

        from ampesm.embeddings import masked_mean_pool

        dim = 4
        # Lote de dos: la primera con 3 residuos reales, la segunda con 5.
        # Tokens: 0 = cls, 2 = eos, 1 = pad, el resto residuos.
        input_ids = torch.tensor([[0, 5, 6, 7, 2, 1, 1], [0, 5, 6, 7, 8, 9, 2]])
        attention_mask = torch.tensor([[1, 1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 1, 1, 1]])
        hidden = torch.arange(2 * 7 * dim, dtype=torch.float32).reshape(2, 7, dim)

        pooled = masked_mean_pool(hidden, attention_mask, input_ids, [0, 2, 1])

        # Referencia calculada a mano: media de las posiciones 1, 2 y 3 de la
        # primera secuencia (las que no son cls, eos ni pad).
        expected_first = hidden[0, 1:4].mean(dim=0)
        assert torch.allclose(pooled[0], expected_first, atol=1e-5)

        expected_second = hidden[1, 1:6].mean(dim=0)
        assert torch.allclose(pooled[1], expected_second, atol=1e-5)

    def test_special_tokens_excluded(self):
        """<cls> y <eos> no entran en la media.

        ESM-2 se entrena solo con enmascaramiento de tokens, sin objetivo a
        nivel de secuencia, asi que <cls> no tiene significado global
        garantizado: incluirlo solo anade una constante que difumina la senal,
        y mas en los peptidos cortos.
        """
        import torch

        from ampesm.embeddings import masked_mean_pool

        input_ids = torch.tensor([[0, 5, 6, 2]])
        attention_mask = torch.ones_like(input_ids)
        hidden = torch.tensor([[[100.0], [1.0], [3.0], [100.0]]])

        pooled = masked_mean_pool(hidden, attention_mask, input_ids, [0, 2, 1])
        assert torch.allclose(pooled, torch.tensor([[2.0]]))  # media de 1 y 3

    def test_all_special_does_not_divide_by_zero(self):
        """Caso degenerado: si todo quedara enmascarado, no debe dar NaN."""
        import torch

        from ampesm.embeddings import masked_mean_pool

        input_ids = torch.tensor([[0, 2]])
        attention_mask = torch.ones_like(input_ids)
        hidden = torch.tensor([[[1.0], [2.0]]])

        pooled = masked_mean_pool(hidden, attention_mask, input_ids, [0, 2, 1])
        assert torch.isfinite(pooled).all()


class TestModelRegistry:
    def test_dims_and_sizes_cover_same_models(self):
        from ampesm.embeddings import MODEL_SIZES_MB

        assert set(MODEL_DIMS) == set(MODEL_SIZES_MB)

    def test_dimensions_grow_with_model_size(self):
        from ampesm.embeddings import MODEL_SIZES_MB

        pairs = sorted((MODEL_SIZES_MB[m], MODEL_DIMS[m]) for m in MODEL_DIMS)
        dims = [d for _, d in pairs]
        assert dims == sorted(dims)
