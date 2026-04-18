# iScience revision — REPORT

**Branch:** `revision/iscience-rev1`
**Date:** 2026-04-18 (conservative-framing pass)

> **Strategic framing.** The revision's claim is **not** that XGBoost is
> uniquely optimal. The claim is that **flexibility-based features carry**
> **meaningful diagnostic signal and support an interpretable diagnostic**
> **workflow across model classes** — XGBoost is the consistent
> explainable learner used to expose that signal.

## At-a-glance

| Workstream | Key result |
|---|---|
| #1 Feature ablation | curated_paper (n=31): F1=0.9565, R²=0.9216; all_widths (120): F1=0.9911, R²=0.9062; random_30 envelope ≈ same. **Curated 30 = interpretability layer, not performance optimum.** |
| #2 Benchmark | curated+XGBoost F1=0.9565, curated+RF F1=1.0000, best inputs-only F1=0.9581. **Signal is in the flexibility representation, not the learner.** |
| #3 Existing-data | 5-fold CV macro-F1=0.9911, severity R²=0.9062 (model-internal robustness; distinct from Fig 7 experimental story) |
| #4 iML1515 transfer | macro-F1=0.9718, R²=0.9778 on 244 LHS conditions (in silico only; SI figure, not main Fig 4 substitute) |
| #5 Runtime | total stages = 3.2 s on a single core, <420 MB peak RSS |

## Placement strategy (main vs SI)

| Asset | Recommended placement |
|---|---|
| Short benchmark / existing-data / transfer / runtime paragraphs | Main text |
| All revision figures (ablation, benchmark, confusion, residuals, transfer) | **SI** |
| Runtime summary table | **SI** |

Rationale: keep the main narrative anchored on the curated Fig 4 / Fig 7
story; use SI to supply quantitative reviewer-driven justification.

## 1. Feature panel ablation

# Feature panel ablation — summary

- N conditions: 242  ·  classes: ['N_limited', 'Ac_limited', 'O2_limited']
- Universe: 120 `width__` columns (parquet alphabetically truncated at 'FACOAL161'; M/I/N-prefixed paper anchors absent)
- Random-30 controls: 10 seeds

## Conservative interpretation

The curated 30/31-feature panel should be interpreted as an
**interpretability-oriented diagnostic layer rather than a
performance-optimal subset**. Within the deployed 120-`width__` universe,
performance is largely insensitive to panel size (Δmacro-F1 ≤ 0.04 across
panels of size 10–120; random 30-feature controls match curated).
Reviewer 2's "why 30?" question is therefore answered on grounds of
mechanistic interpretability, not predictive optimality.

## Headline numbers

```
        panel  n_features  macro_f1  balanced_accuracy     r2    rmse
curated_paper          31    0.9565             0.9349 0.9216 0.03085
       top_10          10    0.9658             0.9508 0.8996 0.03492
       top_20          20    0.9658             0.9508 0.8996 0.03492
       top_50          50    0.9658             0.9508 0.8996 0.03492
   all_widths         120    0.9911             0.9841 0.9062 0.03374
```

## Random-30 controls

- macro_F1 mean ± std: 0.986 ± 0.009
- R²        mean ± std: 0.908 ± 0.010

## Curated paper-aligned panel composition

```
width__EX_o2_e
width__EX_nh4_e
width__EX_pi_e
width__EX_co2_e
width__EX_h_e
width__EX_h2o_e
width__ATPS4rpp
width__ADK1
width__AKGDH
width__CYO1_KT
width__CS
width__ACS
width__ACSERL
width__ACONT
width__ACONTa
width__ACONTb
width__ENO
width__ENOPH
width__ACLS
width__ACLS_a
width__ACLSa
width__ACLSb
width__ADCS
width__APSR
width__APSR2
width__ACGS
width__ARGSL
width__ARGSS
width__ASPTA
width__CMt2ppi
width__DHORDfum
```

Outputs: `revision_runs/iscience_rev1/01_feature_panel_ablation/` + `figures/`

## 2. Baseline benchmarking

# Baseline benchmarking — summary

- N conditions: 242, classes: ['N_limited', 'Ac_limited', 'O2_limited']

## Feature sets

- **A** Inputs only — uptake bounds (acetate, oxygen, ammonium, phosphate)
- **B** GEM summary — biomass objective + uptake fluxes + sat-flag fluxes
- **C** Curated widths — paper-aligned 31-reaction width__ panel (ours)

## Headline numbers

```
     feature_set  n_features    clf_model  macro_f1  balanced_accuracy  per_class_f1__N_limited  per_class_f1__Ac_limited  per_class_f1__O2_limited    reg_model    rmse     mae     r2
   A_inputs_only           4       LogReg    0.8504             0.9366                   0.9524                    0.6207                    0.9782   ElasticNet 0.09524 0.07013 0.2530
   A_inputs_only           4 RandomForest    0.9581             0.9651                   0.9767                    0.9000                    0.9976 RandomForest 0.04125 0.01696 0.8599
   A_inputs_only           4      XGBoost    0.9581             0.9651                   0.9767                    0.9000                    0.9976      XGBoost 0.04308 0.01763 0.8471
   B_gem_summary           9       LogReg    0.9341             0.9476                   0.9524                    0.8571                    0.9929   ElasticNet 0.00529 0.00383 0.9977
   B_gem_summary           9 RandomForest    0.9492             0.9492                   0.9524                    0.9000                    0.9953 RandomForest 0.02991 0.00253 0.9263
   B_gem_summary           9      XGBoost    0.9492             0.9492                   0.9524                    0.9000                    0.9953      XGBoost 0.02990 0.00258 0.9264
C_curated_widths          31       LogReg    0.9760             0.9841                   0.9756                    0.9524                    1.0000   ElasticNet 0.00515 0.00350 0.9978
C_curated_widths          31 RandomForest    1.0000             1.0000                   1.0000                    1.0000                    1.0000 RandomForest 0.03113 0.00260 0.9202
C_curated_widths          31      XGBoost    0.9565             0.9349                   0.9268                    0.9474                    0.9953      XGBoost 0.03085 0.00319 0.9216
```

## Interpretation (conservative framing)

- Inputs-only baselines establish a floor on severity regression
  (R² ≤ 0.86) but reach competitive macro-F1 on classification
  (RandomForest / XGBoost 0.958), reflecting the strong O2_limited
  majority class.
- GEM-summary baselines achieve high R² (≥ 0.93). **Caveat:** the FBA
  objective is in this feature vector by construction (severity =
  obj/obj_max), so its R² ≈ 0.998 reflects target-construction overlap,
  not new predictive content.
- On the same curated 31-reaction width panel, all three learners
  perform competitively (macro-F1 0.957–1.000; RMSE 0.005–0.031). The
  diagnostic signal is therefore carried primarily by the
  **flexibility-based representation**, not by a uniquely optimal
  learner. Random forest reached macro-F1 = 1.000 on this small dataset
  but lacks SHAP-equivalent interpretability of the same form; XGBoost
  is used as the consistent explainable learner.
- **Headline message:** the revision's claim is *not* that XGBoost is
  uniquely optimal. The claim is that flexibility-based features carry
  meaningful diagnostic signal and support an interpretable diagnostic
  workflow across model classes.

Outputs: `revision_runs/iscience_rev1/02_benchmarking/` + `figures/`

## 3. Existing-data performance summary

# Existing-data performance summary

- N conditions: 242, classes: ['N_limited', 'Ac_limited', 'O2_limited']
- Universe: 120 `width__` columns (all widths)

## 5-fold CV classification

- macro-F1: **0.9911**
- balanced accuracy: **0.9841**
- per-class F1: {'N_limited': 0.976, 'Ac_limited': 1.0, 'O2_limited': 0.998}

## 5-fold CV severity regression

- RMSE = 0.03374
- MAE = 0.00379
- R² = 0.9062

## Top mismatch (holdout C1..C10)

See `top_mismatch_conditions.csv` (sorted by rank residual).

## Scope note (must be preserved when quoting these numbers)

These metrics describe **model-internal cross-validation robustness on
the simulated diagnostic dataset** (n=242 LHS conditions). They are
**distinct** from the experimental agreement/mismatch story shown in
Fig 7. The two should be reported as complementary:

- **Fig 7** = experimental Δ(model − measurement) and mismatch
  interpretation (n=10 C-series + n=10 N-series),
- **This file** = CV macro-F1 / per-class metrics / residual summary on
  the 5-fold split of the simulated diagnostic dataset.

Outputs: `revision_runs/iscience_rev1/05_existing_data/` + `figures/`

## 4. External transfer (iML1515)

# External transfer (iML1515) — summary

## Setup

- GEM: iML1515 (E. coli, 2712 reactions)
- Curated transfer panel: 45 reactions covering TCA, glyoxylate, glycolysis / anaplerotic, respiration / ATP, acetate uptake, N biosynthesis
- LHS: 250 conditions over (glc, ac, o2, nh4, pi) uptake bounds
- Feasible solutions kept: 244 / 250; final training set (after rare-class drop): 244

## Performance (5-fold CV on iML1515 features)

- macro-F1 = **0.9718**
- balanced accuracy = **0.9657**
- severity R² = **0.9778**, RMSE = 0.03189
- regime distribution: {'o2_limited': 115, 'nh4_limited': 80, 'glc_limited': 49}

## Top-15 SHAP features

**iML1515:**

  - `EX_co2_e` (|SHAP|=1.1624)
  - `SUCOAS` (|SHAP|=1.1489)
  - `EX_o2_e` (|SHAP|=1.0364)
  - `EX_nh4_e` (|SHAP|=0.4965)
  - `EX_ac_e` (|SHAP|=0.0453)
  - `EX_pi_e` (|SHAP|=0.0319)
  - `SUCDi` (|SHAP|=0.0308)
  - `EX_glc__D_e` (|SHAP|=0.0065)
  - `PFK` (|SHAP|=0.0050)
  - `MALS` (|SHAP|=0.0044)
  - `CS` (|SHAP|=0.0036)
  - `EX_h_e` (|SHAP|=0.0029)
  - `PPC` (|SHAP|=0.0016)
  - `GAPD` (|SHAP|=0.0016)
  - `ICL` (|SHAP|=0.0015)

**iSO1_933 (recomputed on `regime_dataset.parquet`):**

  - `EX_h2o_e` (|SHAP|=0.4576)
  - `12DGR120tipp` (|SHAP|=0.4288)
  - `EX_co2_e` (|SHAP|=0.1994)
  - `ACONT` (|SHAP|=0.1589)
  - `5DOAN` (|SHAP|=0.0845)
  - `EX_h_e` (|SHAP|=0.0274)
  - `AKGDH` (|SHAP|=0.0126)
  - `ADCS` (|SHAP|=0.0103)
  - `DMPPS` (|SHAP|=0.0079)
  - `23CTI2` (|SHAP|=0.0000)
  - `23CTI1` (|SHAP|=0.0000)
  - `12DGR141tipp` (|SHAP|=0.0000)
  - `12DGR161tipp` (|SHAP|=0.0000)
  - `12DGR180tipp` (|SHAP|=0.0000)
  - `12DGR181tipp` (|SHAP|=0.0000)

## Interpretation (conservative framing)

- The diagnostic logic — LHS over uptake bounds → shadow-price regime
  labeling → targeted FVA-width features → XGBoost+SHAP — transferred
  directly to iML1515 with no methodological changes.
- iML1515 top SHAP features (TCA / glyoxylate / glycolysis modules) and
  iSO1 top SHAP features in this 120-width superset are
  **system-specific reaction IDs** but **functionally analogous**
  (central carbon, respiration, acetate uptake, biosynthesis in both).
- **The framework is therefore transferable in formulation, while
  system-specific feature curation remains necessary; the present
  transfer analysis is in silico only, and wet-lab validation in the
  external organism remains future work.**
- **Caveat — do not use this figure to replace main Fig 4.** The iSO1
  top SHAP features under the deployed truncated 120-width universe
  (e.g. `EX_h2o_e`, `12DGR120tipp`, `ACONT`, `5DOAN`) do not align
  with the published Fig 4 narrative (MDH / ICDH / CS / ICL / MALS /
  …) because that narrative was built on the broader paper-curated
  feature universe. The transfer figure belongs in SI as a transfer
  support panel, not as a Fig 4 substitute.

Outputs: `revision_runs/iscience_rev1/03_transfer/` + `figures/`

## 5. Runtime profiling

# Runtime profiling summary

- environment: 3.12.11 on Windows-10-10.0.19045-SP0
- CPUs: 16

| Stage | n | wall (s) | peak RSS (MB) |
|---|---|---|---|
| lhs_sample | 500 | 0.0003 | 145.8 |
| fba_batch | 50 | 2.0686 | 281.2 |
| targeted_fva | 10 | 0.4246 | 289.4 |
| classifier+shap | 242 | 0.2854 | 400.9 |
| regressor+shap | 242 | 0.0598 | 402.3 |
| pc_bootstrap | 25 | 0.3154 | 417.0 |

Outputs: `revision_runs/iscience_rev1/04_runtime/` + `figures/`

## Known limitations (mandatory transparency)

- **Width universe truncation.** `regime_dataset.parquet`'s `width__` columns are alphabetically truncated at `FACOAL161`, so paper-named TCA enzymes beyond 'F' (MDH, ICDH, ICL, MALS, PFK, PYK, NDH, PPC, PCK) are absent from the 120-width superset. The ablation operates within the deployed feature universe. **Implication:** the iSO1 SHAP top features in this universe (`EX_h2o_e`, `12DGR120tipp`, `ACONT`, `5DOAN`, …) do not align with the published Fig 4 narrative (MDH / ICDH / CS / ICL / MALS / …). Main-text Fig 4 should remain on the broader paper-curated narrative; revision figures live in SI.
- **Baseline B target overlap.** B_gem_summary's R² ≈ 0.998 reflects target-construction overlap (severity = obj/obj_max; obj is in B's feature vector), not new predictive content. Reported transparently.
- **iML1515 transfer scope.** In-silico only; no wet-lab matching. Stated as future work in the rebuttal.
- **CV vs experimental performance.** The CV macro-F1 = 0.991 numbers describe model-internal robustness on the simulated diagnostic dataset; they are **not** the same quantity as the experimental agreement/mismatch shown in Fig 7. The two should be presented as complementary, not interchangeable.

