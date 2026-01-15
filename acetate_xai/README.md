## acetate-xai

Genome-scale metabolic digital twin pipeline (FBA + targeted FVA) to support explainable analysis of growth-limiting regimes under acetate-fed conditions.

### What you get in Step 1 (project skeleton)

- `src/acetate_xai/cli.py`: Typer-based CLI (currently includes `--help` and a `version` command)
- `src/acetate_xai/io.py`: placeholders for SBML load + table save helpers
- `src/acetate_xai/config.py`: YAML/JSON config loading helpers
- `tests/test_smoke.py`: smoke test to ensure CLI import/`--help` works

### Local setup (example)

```bash
cd "Stenotrophomonas-causal AI/acetate_xai"
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
pytest -q
python -m acetate_xai.cli --help
```

### Next steps (coming next)

- Implement `simulate fba` and `simulate fva` CLI commands
- Add `configs/` schema for medium + experiment table integration
- Add a single-condition smoke run that produces `results/*.parquet`

