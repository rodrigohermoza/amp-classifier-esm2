.PHONY: help setup lint test check-env data split baselines embeddings train evaluate interpret all clean clean-data

# P-28: el entorno virtual vive FUERA del proyecto a proposito. El proyecto esta
# en ~/Desktop, iCloud sincroniza esa carpeta, y ante un conflicto iCloud no
# fusiona: conserva las dos copias y renombra una a "nombre 2.py". Un venv son
# ~14 000 archivos que Python localiza por nombre exacto, asi que basta ese
# renombrado para dejarlo inservible sin borrar nada. Llego a 4 785 archivos
# duplicados, entre ellos todos los de pytest y las cuatro copias del .pth del
# instalable editable.
#
# Esa era la causa real de P-13 y P-25, que reaparecieron cinco veces porque se
# trataba el sintoma (reinstalar) y no la causa. Mientras el venv estuviera bajo
# ~/Desktop, volveria a pasar.
# `.venv` dentro del proyecto es un ENLACE a esa ruta, no un directorio: asi
# `uv run` a secas, sin pasar por make, tambien resuelve fuera de iCloud.
export UV_PROJECT_ENVIRONMENT ?= $(HOME)/.venvs/ampesm

PY := uv run

help:
	@echo "Objetivos disponibles:"
	@echo "  setup       Crea el entorno con uv e instala dependencias"
	@echo "  lint        Ejecuta ruff"
	@echo "  test        Ejecuta pytest"
	@echo "  check-env   Comprueba que el paquete se importa (P-25, P-28)"
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
	@test -e .venv || ln -s "$(UV_PROJECT_ENVIRONMENT)" .venv
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
	  echo "Entorno: $$UV_PROJECT_ENVIRONMENT"; \
	  dup=$$(find "$$UV_PROJECT_ENVIRONMENT" -name "* [0-9]" -o -name "* [0-9].*" 2>/dev/null | wc -l | tr -d " "); \
	  if [ "$$dup" -gt 0 ]; then \
	    echo "Causa: $$dup archivos renombrados por sincronizacion en la nube (P-28)."; \
	    echo "El entorno esta dentro de una carpeta sincronizada. Muevelo fuera:"; \
	    echo "  rm -rf \"$$UV_PROJECT_ENVIRONMENT\" && make setup"; \
	  else \
	    echo "Arreglo:"; \
	    echo "  uv sync --reinstall-package ampesm"; \
	    echo "  uv sync --extra dev   # el comando anterior quita las deps de desarrollo"; \
	  fi; \
	  exit 1; }
	@if [ -d .venv ] && [ ! -L .venv ]; then \
	  echo "Aviso (P-28): .venv es un directorio real dentro del proyecto, no un"; \
	  echo "enlace. Funciona hoy, pero iCloud lo corrompera. Reconvertirlo:"; \
	  echo "  rm -rf .venv && ln -s \"$$UV_PROJECT_ENVIRONMENT\" .venv"; \
	fi

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
