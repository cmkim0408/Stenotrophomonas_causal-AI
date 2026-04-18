# Project goal

This repository is being revised for an **iScience major revision**.
The goal is to run computational analyses required for revision, not to rewrite the whole project.

# Hard constraints

- Do not modify manuscript, SI, cover-letter, or release-package files unless explicitly asked.
  Specifically off-limits unless the user says otherwise:
  - `_unpacked_main/`, `_unpacked_si/`, `_extract_docx.py` (transient inspection artifacts)
  - `Stenotrophomonas_causal-AI-results-2026-01-16.{tar.gz,zip}` and other release archives
  - `BaseModel.xml`, `Final_genome.fasta` (frozen artifacts referenced by the published draft)
  - Existing committed scripts under `acetate_xai/scripts/` — extend in `code/revision/` instead
- Work only under:
  - `code/revision/`            — new revision scripts (importable from existing `acetate_xai`)
  - `revision_runs/iscience_rev1/` — all outputs (csv/json/parquet/png/pdf)
  - `docs/revision/`            — plan, report, rebuttal inserts
- Reuse existing pipeline code whenever possible:
  - `acetate_xai/src/acetate_xai/` (xai, fva, regime, io, config, medium)
  - `acetate_xai/scripts/train_xgb_shap_regime.py`, `train_xgb_shap_severity.py`
  - `acetate_xai/scripts/build_regime_dataset.py`, `run_fva_batch*.py`
  - Configs: `acetate_xai/configs/anchors.yaml`, `targets_120.json`, `targets_300.json`
- Prefer additive changes over destructive refactors.
- Keep all runs deterministic: fixed `--seed 42` (random-control panels: seed grid 1..10).
- Save every result in machine-readable format (csv/json/parquet) **and** figure format (png/pdf).
- Write a short markdown summary after each completed analysis.

# Compute environment (verified 2026-04-18)

- Python interpreter: `C:/Users/user/Miniforge3/python.exe` (Python 3.12.2)
- Solver: `optlang.glpk_interface` (default; works for FBA/FVA on `BaseModel.xml`)
- Installed: cobra 0.31.1 (user-site overlay) / 0.29.1 (base), xgboost 3.1.3, sklearn 1.8.0, pandas 2.3.3, scipy 1.17.1, shap 0.49.1
- **Known issues to fix before runs:**
  1. `numpy 2.4` conflicts with `numba` (breaks `import shap`) → pin `numpy<=2.3` or upgrade `numba`
  2. `causallearn` is **not installed** — needed for PC bootstrap → `pip install causal-learn`
  3. xgboost version skew: paper says 2.0.0, env has 3.1.3 → existing model JSON files load via the metadata-compat shim already added in commit `de7ce35`
- E. coli **iML1515** loads via `cobra.io.load_model('iML1515')` (2,712 reactions). Use this for the external-transfer workstream.

# Required analyses (priority order)

1. **Feature panel ablation / expansion** — answers Reviewer 2 "왜 30개?" (mandatory 2 + 3)
2. **Benchmarking against baselines** — answers editor + Reviewer (mandatory 4 + 3)
3. **Existing-data performance summaries** — confusion matrix, residuals, mismatch table (mandatory 3 supplement)
4. **External transfer / generalizability** — iML1515 in silico (mandatory 5 + 6)
5. **Runtime / scalability profiling** — wall-clock + env capture (mandatory 9)

# Deliverables (each must exist before Plan-B is reported done)

- `revision_runs/iscience_rev1/REPORT.md`
- `revision_runs/iscience_rev1/metrics_summary.csv`
- `revision_runs/iscience_rev1/figures/` (png+pdf per workstream)
- `docs/revision/rebuttal_insertions.md`

# Workflow rule

Plan first (`docs/revision/iscience_revision_plan.md`), wait for user approval, then implement under `code/revision/` and emit outputs into `revision_runs/iscience_rev1/`. After each workstream completes, append a 3-bullet summary to `REPORT.md`.
