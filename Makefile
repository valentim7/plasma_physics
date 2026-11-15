.PHONY: setup data test lint reproduce note all clean

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

note:  ## build docs/technical-note.pdf from docs/technical-note/note.typ
	uv run python -c "import typst; typst.compile('docs/technical-note/note.typ', output='docs/technical-note.pdf', root='.')"

all: setup data lint test reproduce note

clean:
	rm -rf .venv .pytest_cache .ruff_cache
