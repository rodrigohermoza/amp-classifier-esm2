.PHONY: help setup lint test data split baselines embeddings train evaluate interpret all clean clean-data

PY := uv run

help:
	@echo "Objetivos disponibles:"
	@echo "  setup       Crea el entorno con uv e instala dependencias"
	@echo "  lint        Ejecuta ruff"
	@echo "  test        Ejecuta pytest"
	@echo "  data        Descarga y limpia las fuentes -> data/processed/peptides.csv"
	@echo "  split       Agrupa por identidad de secuencia y divide train/val/test"
	@echo "  baselines   Entrena composicion y fisicoquimicos (LR, RF, XGBoost)"
	@echo "  embeddings  Calcula embeddings ESM-2 y los cachea"
	@echo "  train       Entrena los clasificadores sobre embeddings"
	@echo "  evaluate    Metricas en test con IC bootstrap + figuras"
	@echo "  interpret   UMAP y analisis de errores"
	@echo "  all         data -> split -> baselines -> embeddings -> train -> evaluate -> interpret"
	@echo "  clean       Borra resultados derivados (deja los datos)"
	@echo "  clean-data  Borra data/interim y data/processed (deja data/raw)"

setup:
	uv sync --extra dev
	$(PY) pre-commit install

lint:
	$(PY) ruff check src tests
	$(PY) ruff format --check src tests

test:
	$(PY) pytest

data:
	$(PY) ampesm data

# Borra las predicciones a proposito: una division nueva las invalida, y
# cruzarlas con los splits nuevos produce el fallo de P-10a.
split:
	rm -rf results/metrics/predictions
	$(PY) ampesm split

baselines:
	$(PY) ampesm baselines

embeddings:
	$(PY) ampesm embeddings

train:
	$(PY) ampesm train

evaluate:
	$(PY) ampesm evaluate

interpret:
	$(PY) ampesm interpret

all: data split baselines embeddings train evaluate interpret

clean:
	rm -rf results/figures/*.png results/metrics/*.json results/metrics/*.csv \
	       results/metrics/predictions results/models/*
	touch results/figures/.gitkeep results/metrics/.gitkeep results/models/.gitkeep

clean-data:
	rm -rf data/interim/* data/processed/*
	touch data/interim/.gitkeep data/processed/.gitkeep
