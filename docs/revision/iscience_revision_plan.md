# iScience Revision — Implementation Plan (Plan A)

**Branch:** `revision/iscience-rev1`
**Date:** 2026-04-18
**Status:** Plan only — awaiting user approval before any code is written.

---

## 0. Inventory snapshot

### 0.1 Already available (immediately reusable)

| Asset | Path | Notes |
|---|---|---|
| LHS dataset (FBA + targeted FVA) | `results/regime_dataset.parquet` | 242 conditions × 408 cols. **120 `width__` features** (= the universe for ablation). Labels: O2_limited 211 / N_limited 21 / Ac_limited 10. |
| Holdout predictions | `results/holdout_predictions.csv` | C1–C10 prereg outputs (Ac_limited × 8, O2_limited × 2). |
| Holdout OD600 | `data/holdout_od_results.csv` | Filled experimental OD values. |
| Sensitivity outputs | `results/sensitivity/` | eps 0.5x / 2x / baseline. |
| GEM (this study) | `BaseModel.xml` | iSO1_933 (2,092 rxns, 1,367 mets). GLPK solves; needs medium application before FBA. |
| Curated 30 panel | `acetate_xai/configs/anchors.yaml` | 37 named anchor groups (paper text says 30 reactions/modules). |
| Expanded panels | `acetate_xai/configs/targets_120.json`, `targets_300.json` | 120 / 300 reaction lists already enumerated. |
| Trained classifier (legacy) | `results/figures_draft/Fig04_intracellular_only_regime_xgb/model.json` | Reusable as a starting weight check, not as the canonical model. |
| Pipeline source | `acetate_xai/src/acetate_xai/{xai,fva,regime,io,config,medium}.py` | Importable as a package after `pip install -e acetate_xai`. |
| Training scripts | `acetate_xai/scripts/train_xgb_shap_regime.py`, `train_xgb_shap_severity.py`, `build_regime_dataset.py` | Will be **imported** from new revision scripts, not modified. |

### 0.2 Missing / blocking

| Item | Workstream | Blocking? | Resolution |
|---|---|---|---|
| `causallearn` package | (not used in revision workstreams 1–5) | No | install only if PC stability is re-run |
| `numpy ≤ 2.3` (or `numba` upgrade) | All workstreams that import `shap` | **Yes** | `pip install "numpy<2.4"` once at start |
| iML1515 SBML cache | #4 External transfer | No | `cobra.io.load_model('iML1515')` already works (auto-downloads to user cache) |
| External-system experimental data | #4 (full transfer) | Partial | Workstream proceeds **in silico only**; experimental matching deferred and noted as future work |

### 0.3 Compute env (verified)

- Python 3.12.2 at `C:/Users/user/Miniforge3/python.exe`
- cobra 0.31.1 (user-site overlay), xgboost 3.1.3, scikit-learn 1.8.0, pandas 2.3.3, scipy 1.17.1, shap 0.49.1
- GLPK solver works on `BaseModel.xml`
- Windows 10, OS-reported only — no GPU dependency

---

## 1. Workstreams

### #1 Feature panel ablation / expansion (priority 1)

**Goal:** quantify how the choice of FVA-width feature panel affects regime classification and severity regression. Answer Reviewer 2 "why 30?".

**Panels to compare:**
| Panel | Source | n features |
|---|---|---|
| `curated_30` | `anchors.yaml` → match against `width__` cols | ~30 |
| `top_10` | top-K by SHAP|mean abs| ranked from a curated_30 baseline | 10 |
| `top_20` | top-20 by SHAP from curated_30 baseline | 20 |
| `expanded_50` | top-50 by SHAP from full 120-`width__` superset | 50 |
| `all_widths` | every `width__` column | 120 |
| `random_30` × seeds {1..10} | random sample from full superset | 30 (×10 seeds) |

**For each panel × each task:**
- Classification (regime): `XGBClassifier`, stratified 5-fold CV → macro-F1, balanced accuracy, per-class F1
- Regression (severity / `objective_value` or normalized growth-potential): `XGBRegressor`, 5-fold CV → RMSE, MAE, R²
- Save: `ablation_metrics.csv`, `ablation_summary.md`, `figures/feature_panel_ablation.{png,pdf}` (panel size on x, metric on y, with random-30 envelope)

**Reuse:**
- `train_xgb_shap_regime.py::_to_intracellular_only_feature_set` for the base filter
- `train_xgb_shap_severity.py` regressor config

**New file:** `code/revision/01_feature_panel_ablation.py`

---

### #2 Benchmarking against baselines (priority 1)

**Goal:** show that flexibility-feature pipeline beats simpler alternatives. Answer editor "critical".

**Baselines:**
| ID | Features | Models |
|---|---|---|
| A. Inputs only | acetate_lb, oxygen_lb, ammonium_lb, phosphate_lb (already in parquet) | LogisticRegression, RandomForest, XGBoost |
| B. GEM summaries | objective_value + shadow prices for 4 limiting exchanges | LogisticRegression, RandomForest |
| C. Same-panel simpler ML | curated_30 widths | LogisticRegression (multinomial), RandomForest, ElasticNet (regression), RandomForestRegressor |
| Ours | curated_30 widths + XGBoost | already trained |

**Outputs:** `benchmark_metrics.csv`, `benchmark_summary.md`, `figures/benchmark_comparison.{png,pdf}` (grouped bar: macro-F1 / R² across baselines).

**Reuse:** uptake-bound columns + objective already exist in `regime_dataset.parquet`. Shadow prices may need extraction — check column list before deciding.

**New file:** `code/revision/02_benchmark_baselines.py`

---

### #3 Existing-data performance summary (priority 1)

**Goal:** make the C1–C10 holdout story rigorous (mandatory 3 supplement).

**Outputs from existing files only — no new simulation:**
- Confusion matrix on 5-fold CV predictions over `regime_dataset.parquet` → `figures/confusion_matrix.{png,pdf}`
- Per-class precision / recall / F1 → `performance_metrics.csv`
- Regression residual summary on severity → `regression_residuals.{png,pdf}`, `residual_summary.csv`
- Top-mismatch conditions table from `holdout_predictions.csv` vs `holdout_od_results.csv` (paper Fig 7 red points) → `top_mismatch_conditions.csv`

**Reuse:** `local_validate_holdout.py` already produces a confusion matrix; promote to revision outputs.

**New file:** `code/revision/05_existing_data_performance_summary.py`

---

### #4 External transfer — iML1515 (priority 2)

**Goal:** show diagnostic logic transfers to a different organism. Answer mandatory 5+6.

**Plan (in silico only for v1):**
1. Load iML1515 (`cobra.io.load_model('iML1515')`).
2. Apply analogous condition sweep: vary `EX_ac_e`, `EX_glc__D_e`, `EX_o2_e`, `EX_nh4_e` over LHS n=500.
3. Same labeling logic (largest positive shadow price).
4. Same feature engineering: targeted FVA widths over a curated panel for E. coli central carbon (TCA, glyoxylate shunt, respiratory chain, acetate uptake, biosynthesis modules).
5. Train XGBoost regime classifier + SHAP, compare top features qualitatively against iSO1_933 results.
6. Note as "in silico transfer demo; experimental matching deferred."

**Outputs:** `transfer_metrics.csv`, `transfer_summary.md`, `figures/external_transfer.{png,pdf}` (side-by-side SHAP top features iSO1 vs iML1515).

**New files:**
- `code/revision/03_transfer_external_system.py`
- `code/revision/transfer_targets_ecoli.json` (curated panel; ≈40 reactions covering TCA, glyoxylate, respiration, acetate uptake, key biosynthesis)

---

### #5 Runtime / scalability profiling (priority 3)

**Goal:** quick deliverable. Answer mandatory 9.

**Stages timed:** LHS sampling → FBA batch → targeted FVA → classifier train+SHAP → regressor train+SHAP → (optional PC bootstrap if `causallearn` is installed).

**Outputs:**
- `runtime_summary.csv` (stage, n, wall_seconds, peak_RSS_MB)
- `environment_summary.txt` (CPU model, RAM, OS, Python ver, key package vers, solver)

**New file:** `code/revision/04_profile_runtime.py`

---

## 2. Folder layout (will be created during Plan B)

```
code/revision/
  __init__.py
  01_feature_panel_ablation.py
  02_benchmark_baselines.py
  03_transfer_external_system.py
  04_profile_runtime.py
  05_existing_data_performance_summary.py
  06_make_revision_figures_and_tables.py
  utils.py                       # shared metric helpers, panel resolvers
  transfer_targets_ecoli.json
revision_runs/iscience_rev1/
  01_feature_panel_ablation/
    ablation_metrics.csv
    ablation_summary.md
  02_benchmarking/
    benchmark_metrics.csv
    benchmark_summary.md
  03_transfer/
    transfer_metrics.csv
    transfer_summary.md
  04_runtime/
    runtime_summary.csv
    environment_summary.txt
  05_existing_data/
    performance_metrics.csv
    residual_summary.csv
    top_mismatch_conditions.csv
  figures/
    feature_panel_ablation.{png,pdf}
    benchmark_comparison.{png,pdf}
    external_transfer.{png,pdf}
    confusion_matrix.{png,pdf}
    regression_residuals.{png,pdf}
  REPORT.md
  metrics_summary.csv
docs/revision/
  iscience_revision_plan.md   (this file)
  rebuttal_insertions.md      (final synthesis, written last)
```

---

## 3. Categorisation per the user's A/B/C scheme

**(A) Immediately runnable (no new external data)**
- #1 Feature ablation
- #2 Benchmarking
- #3 Existing-data summary
- #5 Runtime profiling

**(B) Needs minor code additions only**
- #4 External transfer (iML1515): needs `transfer_targets_ecoli.json` + a transfer driver. iML1515 itself auto-downloads.

**(C) Blocked by missing external data**
- Real wet-lab E. coli acetate validation — out of scope for v1; will be flagged in `REPORT.md` as future work.

---

## 4. Step-B execution order (when approved)

1. One-time env fix: `pip install "numpy<2.4" causal-learn` (numba/shap and PC bootstrap)
2. Smoke test: load `regime_dataset.parquet`, train one XGBoost on curated_30, check end-to-end works
3. Implement & run #5 runtime first (cheap; gives env baseline)
4. Implement & run #1 ablation
5. Implement & run #2 benchmarking
6. Implement & run #3 existing-data summary
7. Implement & run #4 external transfer
8. `06_make_revision_figures_and_tables.py` consolidates `metrics_summary.csv` and `REPORT.md`
9. Last pass: `docs/revision/rebuttal_insertions.md` — Results / Methods inserts + point-by-point response, conservatively phrased

---

## 5. Risks & mitigations

| Risk | Mitigation |
|---|---|
| numpy/numba conflict blocks SHAP | Pin `numpy<2.4` at env step 1 |
| Curated 30 vs anchors.yaml count mismatch (37 vs 30) | Resolve by intersecting `anchors.yaml` keywords with `width__` columns; document the 30 actually used |
| Tiny minority classes (Ac=10, N=21) bias CV | Use stratified CV with `n_splits=5`; report per-class metrics; consider class_weight balanced for sklearn baselines |
| iML1515 transfer not 1:1 with iSO1 panel | Build a separate curated E. coli panel; report **functional analogue**, not identical reaction list |
| Long runtime on dense FVA in #4 | Use targeted FVA on the curated E. coli panel only; sample n=500 (smaller than n=2000) |

---

## 6. What I will NOT do without explicit approval

- Touch any file in `_unpacked_main/`, `_unpacked_si/`, manuscript, SI, cover letter, release archives, or `BaseModel.xml`
- Modify existing `acetate_xai/scripts/*.py`
- Push to remote / open PR
- Run anything that overwrites `results/holdout_predictions.csv`, `results/sensitivity/*`, or any file already committed

End of plan.
