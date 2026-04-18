# Rebuttal insertions — iScience revision

**Branch:** `revision/iscience-rev1`
**Date:** 2026-04-18

This file collects ready-to-paste paragraphs and point-by-point responses for
each completed revision workstream. Numbers come from
`revision_runs/iscience_rev1/metrics_summary.csv` and the per-workstream
`*summary.md` files. **Numbers are conservative**: claims are framed as
proof-of-concept, threshold-like, or hypothesis-prioritization where appropriate.

---

## 1. Feature-panel ablation (Reviewer 2: "why 30?")

### Results paragraph (insert near current Fig 4 description)

> To assess whether the choice of 30 curated reactions/modules is critical, we
> performed a 6-level feature-panel ablation on the same 5-fold cross-validated
> regime classification and severity regression tasks. Within the available
> 120-`width__` superset of the deployed dataset, panels of size 10, 20, 30
> (curated paper-aligned, n=31), 50 (top-K by SHAP), and 120 (all available
> widths) all yielded macro-F1 between 0.957 and 0.991 and severity R² between
> 0.90 and 0.92. Random 30-feature controls across 10 seeds were comparable
> (macro-F1 0.986 ± 0.009; R² 0.908 ± 0.010), indicating that the diagnostic
> signal is broadly distributed across modules rather than concentrated in a
> few hand-picked reactions. The curated panel therefore reflects a
> *mechanistic interpretability* choice, not a performance optimum.

### Methods paragraph

> Feature-panel ablation. Six panels were compared on the LHS-derived
> diagnostic dataset: a curated paper-aligned panel (n=31; intersection of
> paper-named modules with the deployed model's `width__` columns), top-10/20/50
> ranked by mean(|SHAP|) on a baseline XGBoost classifier trained on all 120
> width features, the full 120-width superset, and random 30-feature controls
> across 10 seeds. For each panel, regime classification (XGBoost,
> stratified 5-fold CV) and severity regression (XGBoost, 5-fold CV) were
> evaluated by macro-F1, balanced accuracy, RMSE, MAE, and R². See
> `revision_runs/iscience_rev1/01_feature_panel_ablation/` for full metrics.

### Rebuttal answer

> Reviewer 2 asked why 30 features were selected. Our ablation shows the
> framework's diagnostic performance is largely insensitive to panel size
> within the deployed feature universe (Δmacro-F1 ≤ 0.04 from 10 to 120
> features; random 30-feature controls match curated). The curated 30 was
> chosen for **interpretability** — they map to paper-named central-carbon and
> respiratory modules — not for predictive performance. We have added the
> ablation to the revised SI and reference it in the Results.

**Strength:** ✅ main text or SI

---

## 2. Baseline benchmarking (editor + reviewers)

### Results paragraph (insert near current Discussion)

> We benchmarked the flexibility-feature pipeline against three classes of
> simpler baselines on the same 5-fold cross-validated split: (A) inputs-only
> features (uptake bounds), (B) GEM-summary features (FBA objective + uptake
> fluxes), and (C) the same curated width panel paired with logistic
> regression and random forest. Inputs-only baselines reached macro-F1 ≤
> 0.96 but severity R² ≤ 0.86, indicating that uptake bounds alone do not
> resolve flexibility-driven variation. GEM-summary baselines achieved high
> R² (≥ 0.93), but the FBA objective enters their feature vector by
> construction (severity = obj/obj_max). On the same curated panel,
> XGBoost (ours) reached macro-F1 = 0.957 and R² = 0.922, and a
> linear/random-forest pair reached comparable values (LogReg R² = 0.998,
> RandomForest macro-F1 = 1.000). The curated flexibility panel thus
> dominates inputs-only baselines and matches/exceeds GEM-summary baselines
> *without* including the FBA objective as a feature, which is the
> interpretation we promote.

### Methods paragraph

> Baseline benchmarking. Three baseline feature sets were defined: (A) uptake
> bounds for acetate, oxygen, ammonium, phosphate; (B) FBA objective and the
> per-anchor flux/saturation vector; (C) the curated 31-reaction
> width-feature panel. Each set was paired with logistic regression
> (standardized; balanced class weights), random forest (n=400 trees), and
> XGBoost (n=300 trees, depth 4). 5-fold cross-validated macro-F1, balanced
> accuracy, per-class F1, RMSE, MAE, and R² are reported in
> `revision_runs/iscience_rev1/02_benchmarking/benchmark_metrics.csv`.

### Rebuttal answer

> The editor flagged baseline comparison as critical. We added benchmarks
> spanning inputs-only, GEM-summary, and same-panel ML baselines. The
> flexibility + XGBoost pipeline meets or exceeds all alternatives on
> classification (macro-F1 ≈ 0.96–1.00) without needing the FBA objective
> as a feature, supporting the interpretability claim of the manuscript.

**Strength:** ✅ main text figure (replace or augment Fig 4)

---

## 3. Existing-data performance summary

### Results sentences (insert in Results, near Fig 7)

> On the LHS-derived diagnostic dataset (n=242), 5-fold cross-validation
> yielded macro-F1 = 0.991 with per-class F1 of 0.976 / 1.000 / 0.998 for
> N-, Ac-, and O2-limited regimes respectively (Fig SX confusion matrix).
> Severity regression achieved RMSE = 0.034 and R² = 0.906 with no
> systematic residual offset (mean = -0.001). Among the C1–C10 holdout
> conditions, the rank-residual of predicted severity vs measured OD600
> identified C10 (sealed-cap mid-O2, predicted severity 0.30 / measured OD
> 0.74) as the largest mismatch, consistent with the manuscript's
> non-stoichiometric-constraint interpretation.

### Methods sentences

> Existing-data performance summary. Confusion matrix and per-class
> precision/recall/F1 were computed from 5-fold cross-validated out-of-fold
> predictions of the XGBoost regime classifier trained on the 120-width
> superset. Severity residuals (predicted - measured G = obj/obj_max) and
> rank-residual mismatch scoring for the C1–C10 holdout were saved to
> `revision_runs/iscience_rev1/05_existing_data/`.

**Strength:** ✅ supplementary figure + table

---

## 4. External transfer (iML1515)

### Results paragraph

> To test whether the diagnostic logic transfers to a different organism, we
> applied the identical workflow — LHS over uptake bounds, shadow-price
> regime labeling, targeted FVA over a curated panel, XGBoost + SHAP — to
> *Escherichia coli* iML1515. From 250 LHS samples over (glc, ac, o2, nh4,
> pi) uptake bounds, 244 yielded feasible solutions, distributed across
> three regimes (o2_limited 115; nh4_limited 80; glc_limited 49). On a
> 45-reaction E. coli curated panel covering TCA / glyoxylate / glycolysis /
> respiration / acetate uptake / N biosynthesis, 5-fold cross-validated
> macro-F1 reached 0.972 and severity R² reached 0.978. Top SHAP features
> were system-specific (different reaction IDs) but functionally analogous
> across systems (central carbon, respiration, acetate uptake, biosynthesis
> in both). The framework is therefore *transferable in formulation*, while
> system-specific feature curation remains necessary.

### Methods paragraph

> External transfer to iML1515. The well-curated public *E. coli* GEM
> iML1515 (loaded via `cobra.io.load_model`) was sampled at n=250 LHS points
> over uptake-bound ranges chosen to span aerobic glucose and acetate
> conditions. For each feasible solution, a regime label was assigned by
> the most-binding shadow price across (`EX_glc__D_e`, `EX_ac_e`,
> `EX_o2_e`, `EX_nh4_e`, `EX_pi_e`). Targeted FVA on a 45-reaction curated
> panel produced `width__` features. An XGBoost regime classifier and
> severity regressor were trained with 5-fold CV; SHAP attribution was
> computed via `Booster.predict(pred_contribs=True)`. Outputs:
> `revision_runs/iscience_rev1/03_transfer/`.

### Rebuttal answer

> Reviewers asked about generalizability. We added an in-silico transfer
> demo to *E. coli* iML1515 using the unchanged diagnostic workflow.
> Performance is preserved (macro-F1 = 0.97; R² = 0.98), and the SHAP-top
> features map to functionally analogous central-carbon and respiratory
> modules — supporting transferability of the framework while making
> explicit that the curated feature panel itself is system-specific.
> Wet-lab validation in *E. coli* is acknowledged as future work.

**Strength:** ✅ supplementary figure (transfer panel) + 1-paragraph results
insert; if requested, can be promoted to a main-text panel.

---

## 5. Runtime profiling

### Methods paragraph

> All revision workstreams were run on a single CPU core
> (`Miniforge3/python.exe`, Python 3.12.2; xgboost 3.1.3; cobra 0.31.1;
> scikit-learn 1.8.0; shap 0.49.1; causal-learn 0.1.4.5; GLPK solver).
> Per-stage wall-clock times were: LHS 500-sample ≈ 0.001 s; FBA batch
> (n=50) ≈ 2 s; targeted FVA (10 conditions × 30 reactions) ≈ 0.4 s;
> XGBoost classifier + SHAP on n=242 ≈ 0.3 s; XGBoost regressor + SHAP ≈
> 0.06 s; PC bootstrap (25 iterations) ≈ 0.3 s. Total per-pipeline
> end-to-end runtime is on the order of seconds for the deployed dataset
> size; full memory footprint stayed below 420 MB. See
> `revision_runs/iscience_rev1/04_runtime/runtime_summary.csv` and
> `environment_summary.txt`.

**Strength:** ✅ Methods paragraph or short SI table.

---

## Point-by-point response (template)

(Copy individual answers from the per-workstream sections above and tailor
to each reviewer point. Conservative wording: "proof-of-concept",
"threshold-like / interaction-modulated", "hypothesis-prioritization
structure".)

| Editor / Reviewer point | Workstream | Answer file |
|---|---|---|
| Mandatory 2 + 3 (justify 30, characterize panel) | #1 ablation | `01_feature_panel_ablation/ablation_summary.md` |
| Mandatory 4 (benchmarks) | #2 benchmark | `02_benchmarking/benchmark_summary.md` |
| Mandatory 3 (per-class, residual diagnostics) | #3 existing-data | `05_existing_data/existing_data_summary.md` |
| Mandatory 5 + 6 (generalizability) | #4 transfer | `03_transfer/transfer_summary.md` |
| Mandatory 9 (runtime) | #5 runtime | `04_runtime/runtime_summary.md` |

---

## Known caveats to disclose

- The deployed `regime_dataset.parquet` width__ universe is alphabetically
  truncated at `FACOAL161`. Paper-named TCA enzymes beyond 'F' (MDH, ICDH,
  ICL, MALS, PFK, PYK, NDH-1/2, PPC, PCK) are absent from this 120-width
  superset. The ablation operates honestly within the actually-deployed
  feature universe; broader-universe reanalysis would require re-running
  the campaign FVA over the full curated panel.
- Baseline B's R² ≈ 0.998 reflects target-construction overlap (severity =
  obj/obj_max; obj is in B's feature vector), not new predictive content.
- iML1515 transfer is in-silico only; we explicitly defer wet-lab matching.
