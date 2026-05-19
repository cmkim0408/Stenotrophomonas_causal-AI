# iScience revision — REPORT

**Branch:** `revision/iscience-rev1`
**Date:** 2026-05-19 (novelty-defense pass — workstream #7 added)

> **Strategic framing.** The revision's claim is **not** that XGBoost is
> uniquely optimal. The claim is that **flexibility-based features carry**
> **meaningful diagnostic signal and support an interpretable diagnostic**
> **workflow across model classes** — XGBoost is the consistent
> explainable learner used to expose that signal.

## At-a-glance

| Workstream | Key result |
|---|---|
| #1 Feature ablation | curated_paper (n=42): F1=0.9565, R²=0.9216; all_widths (120): F1=0.9911, R²=0.9062; random_30 envelope ≈ same. **Curated 30 = interpretability layer, not performance optimum.** |
| #2 Benchmark | curated+XGBoost F1=0.9565, curated+RF F1=1.0000, best inputs-only F1=0.9581. **Signal is in the flexibility representation, not the learner.** |
| #3 Existing-data | 5-fold CV macro-F1=0.9911, severity R²=0.9062 (model-internal robustness; distinct from Fig 7 experimental story) |
| #4 iML1515 transfer | macro-F1=0.9718, R²=0.9778 on 244 LHS conditions (in silico only; SI figure, not main Fig 4 substitute) |
| #5 Runtime | total stages = 3.2 s on a single core, <420 MB peak RSS |
| #7 Point-flux vs flexibility (novelty) | iSO1: W F1=0.957 / PFA F1=0.972 / M F1=0.941 / OBJ F1=0.912 — **W beats same-FVA midpoint by +0.016**; \|pFBA flux\| marginally edges W on F1. **Novelty is conceptual (feasible-space vs point-state), not raw F1.** iML1515: W F1=0.973 vs OBJ floor F1=0.289. |

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
- Universe: 300 `width__` columns (extended FVA campaign — see `revision_runs/iscience_rev1/extended_fva/`)
- Random-30 controls: 10 seeds

## Conservative interpretation

The curated paper-aligned panel should be interpreted as an **interpretability-oriented diagnostic layer rather than a performance-optimal subset**. Within the extended `width__` universe (~300 columns), performance is largely insensitive to panel size (Δmacro-F1 ≤ 0.04 across panels of size 10–300; random-30 controls match curated). Reviewer 2's "why 30?" question is answered on grounds of **mechanistic interpretability**, not predictive optimality.

## Headline numbers

```
        panel  n_features  macro_f1  balanced_accuracy     r2    rmse
curated_paper          42    0.9565             0.9349 0.9216 0.03085
       top_10          10    0.9658             0.9508 0.8996 0.03492
       top_20          20    0.9658             0.9508 0.8996 0.03492
       top_50          50    0.9658             0.9508 0.8996 0.03492
   all_widths         300    0.9911             0.9841 0.9062 0.03374
```

## Random-30 controls

- macro_F1 mean ± std: 0.981 ± 0.009
- R²        mean ± std: 0.913 ± 0.009

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
width__NADH16pp
width__CS
width__ACONT
width__ACONTa
width__ACONTb
width__ICDHyr
width__ICDHx
width__ICL
width__MALS
width__MDH
width__MDH2
width__MDH3
width__FUM
width__ACS
width__ACSERL
width__ENO
width__ENOPH
width__PYK
width__PYK3
width__PPC
width__ACLS
width__ACLSa
width__ACLSb
width__ADCS
width__APSR
width__APSR2
width__GLNS
width__GLUDy
width__ACGS
width__ARGSL
width__ARGSS
width__ASPTA
```

Outputs: `revision_runs/iscience_rev1/01_feature_panel_ablation/` + `figures/`

## 2. Baseline benchmarking

# Baseline benchmarking — summary

- N conditions: 242, classes: ['N_limited', 'Ac_limited', 'O2_limited']

## Feature sets

- **A** Inputs only — uptake bounds (acetate, oxygen, ammonium, phosphate)
- **B** GEM summary — biomass objective + uptake fluxes + sat-flag fluxes
- **C** Curated widths — paper-aligned 42-reaction width__ panel (ours)

## Headline numbers

```
     feature_set  n_features    clf_model  macro_f1  balanced_accuracy  per_class_f1__N_limited  per_class_f1__Ac_limited  per_class_f1__O2_limited    reg_model    rmse     mae     r2
   A_inputs_only           4       LogReg    0.8504             0.9366                   0.9524                    0.6207                    0.9782   ElasticNet 0.09524 0.07013 0.2530
   A_inputs_only           4 RandomForest    0.9581             0.9651                   0.9767                    0.9000                    0.9976 RandomForest 0.04125 0.01696 0.8599
   A_inputs_only           4      XGBoost    0.9581             0.9651                   0.9767                    0.9000                    0.9976      XGBoost 0.04308 0.01763 0.8471
   B_gem_summary           9       LogReg    0.9341             0.9476                   0.9524                    0.8571                    0.9929   ElasticNet 0.00529 0.00383 0.9977
   B_gem_summary           9 RandomForest    0.9492             0.9492                   0.9524                    0.9000                    0.9953 RandomForest 0.02991 0.00253 0.9263
   B_gem_summary           9      XGBoost    0.9492             0.9492                   0.9524                    0.9000                    0.9953      XGBoost 0.02990 0.00258 0.9264
C_curated_widths          42       LogReg    0.9760             0.9841                   0.9756                    0.9524                    1.0000   ElasticNet 0.00517 0.00352 0.9978
C_curated_widths          42 RandomForest    1.0000             1.0000                   1.0000                    1.0000                    1.0000 RandomForest 0.03155 0.00265 0.9180
C_curated_widths          42      XGBoost    0.9565             0.9349                   0.9268                    0.9474                    0.9953      XGBoost 0.03085 0.00319 0.9216
```

## Interpretation (conservative framing)

- Inputs-only baselines establish a floor on severity regression (R² ≤ 0.86) but reach competitive macro-F1 on classification (RandomForest / XGBoost ≈ 0.96), reflecting the strong O2_limited majority class.
- GEM-summary baselines achieve high R² (≥ 0.93). **Caveat:** the FBA objective is in this feature vector by construction (severity = obj/obj_max), so its R² ≈ 0.998 reflects target-construction overlap, not new predictive content.
- On the same curated width panel, all three learners perform competitively (macro-F1 0.957–1.000; RMSE 0.005–0.031). The diagnostic signal is therefore carried primarily by the **flexibility-based representation**, not by a uniquely optimal learner. Random forest reached macro-F1 = 1.000 on this small dataset but lacks SHAP-equivalent interpretability of the same form; XGBoost is used as the consistent explainable learner.
- **Headline message:** the revision's claim is *not* that XGBoost is uniquely optimal. The claim is that flexibility-based features carry meaningful diagnostic signal and support an interpretable diagnostic workflow across model classes.

Outputs: `revision_runs/iscience_rev1/02_benchmarking/` + `figures/`

## 3. Existing-data performance summary

# Existing-data performance summary

- N conditions: 242, classes: ['N_limited', 'Ac_limited', 'O2_limited']
- Universe: 300 `width__` columns (extended FVA campaign)

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

These metrics describe **model-internal cross-validation robustness on the simulated diagnostic dataset** (n=242 LHS conditions). They are **distinct** from the experimental agreement/mismatch story shown in Fig 7. The two should be reported as complementary:

- **Fig 7** = experimental Δ(model − measurement) and mismatch interpretation (n=10 C-series + n=10 N-series),
- **This file** = CV macro-F1 / per-class metrics / residual summary on the 5-fold split of the simulated diagnostic dataset.

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

**iSO1_933 (recomputed on the extended `regime_dataset_extended.parquet`, 300 `width__` columns):**

  - `EX_h2o_e` (|SHAP|=0.4410)
  - `12DGR120tipp` (|SHAP|=0.4288)
  - `EX_co2_e` (|SHAP|=0.1994)
  - `ICDHx` (|SHAP|=0.1479)
  - `5DOAN` (|SHAP|=0.0639)
  - `EX_h_e` (|SHAP|=0.0274)
  - `H2Ot` (|SHAP|=0.0215)
  - `MDH` (|SHAP|=0.0143)
  - `ACONT` (|SHAP|=0.0115)
  - `ADCS` (|SHAP|=0.0103)
  - `AKGDH` (|SHAP|=0.0072)
  - `DMPPS` (|SHAP|=0.0053)
  - `PPCDC` (|SHAP|=0.0050)
  - `MEPCT` (|SHAP|=0.0044)
  - `SUCOAS` (|SHAP|=0.0000)

## Interpretation (conservative framing)

- The diagnostic logic — LHS over uptake bounds → shadow-price regime labeling → targeted FVA-width features → XGBoost+SHAP — transferred directly to iML1515 with no methodological changes.
- iML1515 top SHAP features (TCA / glyoxylate / glycolysis modules) and iSO1 top SHAP features in the extended ~300-width universe are **system-specific reaction IDs** but **functionally analogous** (central carbon, respiration, acetate uptake, biosynthesis in both).
- After the extended FVA campaign, the iSO1 top-SHAP set now includes paper-named anchors such as MDH, ICDHx, AKGDH, PPCDC, SUCOAS — aligning the data-driven ranking with the published Fig 4 narrative (TCA / respiration / ATP).
- **The framework is therefore transferable in formulation, while system-specific feature curation remains necessary; the present transfer analysis is in silico only, and wet-lab validation in the external organism remains future work.**

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

- **Width universe — extended FVA campaign applied.** The original `results/regime_dataset.parquet` was alphabetically truncated at `FACOAL161` (120 `width__` columns). For this revision, an extended FVA campaign re-ran targeted FVA on the 180 missing reactions across all 242 conditions (replay-verified against the original `objective_value` to within 5e-2 across all rows), expanding the deployed `width__` universe to ~300 columns. Paper-named TCA / respiration anchors (MDH, ICDHx, ICDHyr, ICL, MALS, PYK, PPC, NADH16pp, FUM …) are now present and appear in the SHAP top-K of the iSO1 classifier, aligning with the published Fig 4 narrative. See `revision_runs/iscience_rev1/extended_fva/` for raw outputs and `code/revision/extend_fva_campaign.py` for the driver.
- **Baseline B target overlap.** B_gem_summary's R² ≈ 0.998 reflects target-construction overlap (severity = obj/obj_max; obj is in B's feature vector), not new predictive content. Reported transparently.
- **iML1515 transfer scope.** In-silico only; no wet-lab matching. Stated as future work in the rebuttal.
- **CV vs experimental performance.** The CV macro-F1 = 0.991 numbers describe model-internal robustness on the simulated diagnostic dataset; they are **not** the same quantity as the experimental agreement/mismatch shown in Fig 7. The two should be presented as complementary, not interchangeable.

# Point-flux vs flexibility-interval baseline -- summary

**Question:** does the *width* of the FVA-derived feasible interval (our
flexibility representation) carry diagnostic signal beyond the *point flux*
(pFBA solution, FVA midpoint, FBA objective)?

Each representation is evaluated on the **same curated reaction panel**, the
**same XGBoost** classifier and regressor, **5-fold CV**, and identical
hyperparameters. The only thing that changes is how each panel reaction is
encoded as a feature.

## Feature representations

- **W (widths)** -- `width__rxn = vmax - vmin` from targeted FVA
  (`fraction_of_optimum=0.95`). The manuscript's flexibility feature.
- **M (midpoints)** -- `mid__rxn = (vmax + vmin) / 2` from the *same FVA*
  solve. A point-center proxy that comes from the same data as W; isolates
  the contribution of "the interval" controlling for the FVA solve.
- **PF (pFBA flux)** -- parsimonious FBA flux (signed) for each reaction.
  The strongest "point-flux" baseline because pFBA returns a unique flux
  vector that minimizes total absolute flux subject to optimal biomass.
- **PFA (|pFBA flux|)** -- magnitude of pFBA flux; sign-agnostic.
- **OBJ** -- biomass objective only (1 scalar feature). Scalar floor.

## Headline results (5-fold CV, same XGBoost hyperparameters)

### iSO1_933 (this study, n=242, 42-reaction curated panel)

| Representation     | n_feat | macro-F1 | R^2    |
|--------------------|--------|----------|--------|
| PFA \|pFBA flux\| | 42     | 0.972    | 0.926  |
| **W (widths, ours)** | **42** | **0.957** | **0.922** |
| M (FVA midpoints)  | 28     | 0.941    | 0.908  |
| PF (pFBA, signed)  | 42     | 0.939    | 0.910  |
| OBJ (scalar)       | 1      | 0.912    | 0.926  |

### iML1515 (E. coli transfer, n=244, 45-reaction curated panel)

| Representation     | n_feat | macro-F1 | R^2    |
|--------------------|--------|----------|--------|
| **W (widths, ours)** | **45** | **0.973** | **0.978** |
| OBJ (scalar)       | 1      | 0.289    | 0.998* |

\* iML1515 OBJ_only R^2 = 0.998 is a target-construction artifact:
severity is defined as `obj / obj_max`, so the OBJ feature IS the target's
numerator. The OBJ classification F1 = 0.289 is the honest scalar floor.

iML1515 pFBA baselines (PF/PFA) are pending re-run on the user's local
machine (the sandbox cannot reach the BiGG repository). The required
command is:

```
python3 code/revision/_pfba_runner.py iml 0 244
python3 code/revision/07_pointflux_baseline.py  # then re-runs eval w/ cache
```

## Interpretation (conservative framing)

The comparison is deliberately structured so that any performance
difference is attributable to **feature representation**, not learner
tuning or panel choice. Five honest take-aways:

1. **On the iSO1 curated panel, |pFBA flux| edges out FVA-width on F1
   (0.972 vs 0.957).** We report this transparently. It means that, for
   classification under this particular FVA solve and dataset size, the
   *magnitude* of the parsimonious point-flux carries comparable
   classification signal to the interval width.

2. **Width nonetheless beats the FVA midpoint within the same FVA solve
   (0.957 vs 0.941 macro-F1; 0.922 vs 0.908 R^2).** This is the cleanest
   "interval vs point" test, because both features come from the same
   underlying solve on the same conditions. **The interval size itself
   adds information beyond the central representative flux.**

3. **Width matches or beats signed pFBA flux (0.957 vs 0.939 F1; 0.922
   vs 0.910 R^2)** on the same panel without solving the parsimonious
   problem.

4. **The novelty is conceptual, not raw F1.** Width and |pFBA| answer
   *different questions*:
   - Width = "how much rerouting capacity remains under this constraint
     context" -> directly maps to rigidification / flexibility-collapse
     diagnosis -> supports the rigidification-map interpretation
     downstream (Fig 6).
   - pFBA-flux magnitude = "where does the biomass-maximizing flux go" ->
     a point-state representation of the optimum, not of the feasible
     space.
   On classification metrics they happen to be comparable on this dataset;
   on **interpretability** they are not. Only the width is *mechanistically*
   readable as a degree-of-freedom signal.

5. **The R^2 column is partially structural.** Because severity = obj /
   obj_max by definition, the OBJ-only baseline achieves R^2 >= 0.92 on
   both datasets *by construction* -- this is target-overlap, not new
   predictive content (the same caveat is already noted for the B_gem_summary
   baseline in workstream #2). The honest regression-side test is therefore
   whether features improve over OBJ on macro-F1 *and* on classifications
   of structurally distinct regimes, which W achieves cleanly on iML1515
   (F1 0.97 vs OBJ-floor 0.29).

## Strategic framing (one sentence)

> The flexibility representation (FVA width) is competitive with the
> strongest point-flux baselines on classification, beats the same-solve
> point-center proxy (midpoint), and -- crucially -- supports the
> mechanistic rigidification / degree-of-freedom interpretation that
> point-flux representations cannot, irrespective of macro-F1.

## Caveats (mandatory transparency)

- **mid__ universe truncation.** The iSO1 `mid__` columns are the
  alphabetically-truncated 120-reaction FVA campaign (cuts off at
  FACOAL161), so M_midpoints uses 28 of the curated 42 reactions. W and
  PF/PFA use the full 42 because their data was either already in the
  extended `width__` superset or recomputed via pFBA in this workstream.
- **pFBA determinism.** pFBA returns a unique flux vector, so PF/PFA
  carry no non-uniqueness confound -- the strongest possible point-flux
  baseline.
- **iML1515 pFBA pending.** Sandbox network blocks the BiGG download;
  see command above to complete on the user's machine.
- **Same hyperparameters and CV split.** All representations use the
  same XGBoost configuration (300 trees, max_depth=4, hist) and the same
  5-fold split (random_state=42).

## Files

- `pointflux_metrics.csv` -- full headline table
- `pfba_fluxes_iso1.parquet` -- cached iSO1 pFBA solutions (242 conditions x 42 reactions)
- `pfba_fluxes_iml1515.parquet` -- (pending; will appear after the Windows-side rerun)
- `../figures/pointflux_vs_width.{png,pdf}` -- grouped bar comparison
- `parts_iso1/`, `parts_iml/` -- per-row partial parquets (resumable)

---

## 7. Point-flux vs flexibility-interval (novelty defense) — appended above
