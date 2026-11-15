.PHONY: setup data test lint reproduce all clean

setup:  ## create .venv with Python 3.12 and install dependencies
	uv sync

data:  ## check that data/processed/ matches the published checksums
	cd data/processed && shasum -a 256 -c SHA256SUMS

lint:
	uv run ruff check .
	uv run ruff format --check .

test:
	uv run pytest -q

reproduce:  ## write reports/results.md and reports/figures/
	uv run python -m plasma_physics.report

all: setup data lint test reproduce

clean:
	rm -rf .venv .pytest_cache .ruff_cache
