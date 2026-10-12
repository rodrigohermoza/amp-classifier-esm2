.PHONY: help setup lint test check-env data split baselines embeddings train evaluate interpret all clean clean-data

PY := uv run

help:
	@echo "Objetivos disponibles:"
	@echo "  setup       Crea el entorno con uv e instala dependencias"
	@echo "  lint        Ejecuta ruff"
	@echo "  test        Ejecuta pytest"
	@echo "  check-env   Comprueba que el paquete se importa (P-25)"
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

# P-25: el .pth del instalable editable deja de procesarse cada cierto tiempo y
# el ejecutable `ampesm` muere con ModuleNotFoundError antes de tocar nada. Ha
# pasado tres veces, la ultima tumbando el primer `make data` real. El
# pythonpath de pyproject.toml protege a pytest pero NO al ejecutable, asi que
# los tests seguian pasando con la tuberia muerta.
#
# La comprobacion vive aqui, delante de cada etapa, y no en `ampesm info`:
# cuando el import falla, `ampesm info` tampoco arranca. Una guarda dentro del
# programa que no puede ejecutarse no es una guarda (la leccion de P-20).
check-env:
	@$(PY) python -c "import ampesm" 2>/dev/null || { \
	  echo "Error: el paquete 'ampesm' no se importa en este entorno (P-25)."; \
	  echo "Arreglo:"; \
	  echo "  uv sync --reinstall-package ampesm"; \
	  echo "  uv sync --extra dev   # el comando anterior quita las deps de desarrollo"; \
	  exit 1; }

data: check-env
	$(PY) ampesm data

# Borra las predicciones a proposito: una division nueva las invalida, y
# cruzarlas con los splits nuevos produce el fallo de P-10a.
split: check-env
	rm -rf results/metrics/predictions
	$(PY) ampesm split

baselines: check-env
	$(PY) ampesm baselines

embeddings: check-env
	$(PY) ampesm embeddings

train: check-env
	$(PY) ampesm train

evaluate: check-env
	$(PY) ampesm evaluate

interpret: check-env
	$(PY) ampesm interpret

all: data split baselines embeddings train evaluate interpret

clean:
	rm -rf results/figures/*.png results/metrics/*.json results/metrics/*.csv \
	       results/metrics/predictions results/models/*
	touch results/figures/.gitkeep results/metrics/.gitkeep results/models/.gitkeep

clean-data:
	rm -rf data/interim/* data/processed/*
	touch data/interim/.gitkeep data/processed/.gitkeep
