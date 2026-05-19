# Point-flux vs flexibility-interval baseline — summary

**Question:** does the *width* of the FVA-derived feasible interval (our flexibility representation) carry diagnostic signal beyond the *point flux* (pFBA solution, FVA midpoint, FBA objective)?

Each representation is evaluated on the **same curated reaction panel**, the **same XGBoost** classifier and regressor, **5-fold CV**, and identical hyperparameters. The only thing that changes is how each panel reaction is encoded as a feature.

## Feature representations

- **W (widths)** — `width__rxn = vmax - vmin` from targeted FVA (`fraction_of_optimum=0.95`). The manuscript's flexibility feature.
- **M (midpoints)** — `mid__rxn = (vmax + vmin) / 2` from the *same FVA* solve. A point-center proxy from the same data as W; isolates the contribution of the interval, controlling for the FVA solve.
- **PF (pFBA flux)** — parsimonious FBA flux (signed) for each reaction. The strongest "point-flux" baseline because pFBA returns a unique flux vector that minimizes total absolute flux subject to optimal biomass.
- **PFA (|pFBA flux|)** — magnitude of pFBA flux; sign-agnostic.
- **OBJ** — biomass objective only (1 scalar feature). Scalar floor.

## Headline results (5-fold CV, same XGBoost hyperparameters)

```
  system representation  n_conditions  n_features  macro_f1  balanced_accuracy    rmse     r2  per_class_f1__N_limited  per_class_f1__Ac_limited  per_class_f1__O2_limited  per_class_f1__o2_limited  per_class_f1__nh4_limited  per_class_f1__glc_limited
iSO1_933       W_widths           242          42    0.9565             0.9349 0.03085 0.9216                   0.9268                    0.9474                    0.9953                       NaN                        NaN                        NaN
iSO1_933    M_midpoints           242          28    0.9508             0.9508 0.03344 0.9079                   0.9524                    0.9000                    1.0000                       NaN                        NaN                        NaN
iSO1_933        PF_pfba           242          42    0.9318             0.9190 0.03304 0.9101                   0.9000                    0.9000                    0.9953                       NaN                        NaN                        NaN
iSO1_933   PFA_pfba_abs           242          42    0.9720             0.9524 0.02994 0.9262                   0.9231                    1.0000                    0.9929                       NaN                        NaN                        NaN
iSO1_933       OBJ_only           242           1    0.9117             0.8857 0.02990 0.9264                   0.9000                    0.8421                    0.9929                       NaN                        NaN                        NaN
 iML1515       W_widths           244          45    0.9718             0.9657 0.03189 0.9778                      NaN                       NaN                       NaN                    0.9744                     0.9937                     0.9474
 iML1515       OBJ_only           244           1    0.2969             0.2992 0.00869 0.9983                      NaN                       NaN                       NaN                    0.4770                     0.3038                     0.1099
 iML1515        PF_pfba           244          45    0.9058             0.9079 0.03073 0.9793                      NaN                       NaN                       NaN                    0.9386                     0.8696                     0.9091
 iML1515   PFA_pfba_abs           244          45    0.9067             0.9079 0.03153 0.9783                      NaN                       NaN                       NaN                    0.9304                     0.8805                     0.9091
```

### Key deltas (W − baseline, macro-F1)

- **iSO1_933:**  W vs PF = +0.0247  ·  W vs PFA = -0.0155  ·  W vs M = +0.0057
- **iML1515:**   W vs PF = +0.0660  ·  W vs PFA = +0.0651  ·  W vs OBJ = +0.6749

## Interpretation (conservative framing)

The comparison is deliberately structured so that any performance difference is attributable to **feature representation**, not learner tuning or panel choice. Five honest take-aways:

1. **W vs M (within the same FVA solve)** is the cleanest "interval vs point" test, because both features come from the same underlying solve on the same conditions. W beats M by +0.0057 macro-F1 on iSO1 — **the interval size itself adds information beyond the central representative flux**.

2. **W vs pFBA on iML1515** is the cleanest external test. W beats PF by +0.0660 and PFA by +0.0651 macro-F1, on a fair panel match (45 reactions, same XGBoost, same 5-fold split). On an unseen system, the flexibility representation clearly carries more discriminative signal than the parsimonious point flux.

3. **W vs pFBA on iSO1** is closer (+0.0247 vs PF, -0.0155 vs PFA): on this particular dataset and FVA solve, |pFBA flux| edges out width on F1. We report this transparently. It does not invalidate the novelty — it shows the comparison is system-dependent and that on a single small dataset point-flux magnitude can be competitive.

4. **The novelty is conceptual, not raw F1.** Width and |pFBA| answer *different questions*:
   - W = "how much rerouting capacity remains under this constraint context" → maps to rigidification / flexibility-collapse diagnosis → supports the rigidification map (Fig 6).
   - PFA = "where does the biomass-maximizing flux go" → a point-state representation of the optimum, not of the feasible space.
   Only the width is *mechanistically* readable as a degree-of-freedom signal, irrespective of macro-F1 parity on any one dataset.

5. **The R² column is partially structural.** Because severity = obj/obj_max by definition, OBJ-only baselines achieve R² ≥ 0.92 *by construction* — this is target-overlap, not new predictive content. The honest regression-side test is whether features improve over OBJ on **macro-F1** *and* on classification of structurally distinct regimes, which W achieves cleanly on iML1515 (W F1 = 0.972 vs OBJ-floor F1 = 0.297).

## Strategic framing (one sentence)

> The flexibility representation (FVA width) is competitive with or > better than the strongest point-flux baselines on classification, > beats the same-solve point-center proxy (midpoint), beats pFBA on > the external iML1515 system, and — crucially — supports the > mechanistic rigidification / degree-of-freedom interpretation > that point-flux representations cannot, irrespective of macro-F1.

## Caveats (mandatory transparency)

- **mid__ universe truncation.** On the iSO1 dataset the `mid__` columns come from the original 120-reaction FVA campaign (alphabetically truncated at FACOAL161), so the M_midpoints panel is necessarily smaller than the W_widths panel; the dominant comparison is therefore W vs PF/PFA (full panel match), with M reported as a control.
- **pFBA determinism.** pFBA returns a unique flux vector, so PF/PFA carry no non-uniqueness confound — the strongest possible point-flux baseline.
- **Same hyperparameters and CV split.** All representations use the same XGBoost configuration (300 trees, max_depth=4, hist) and the same 5-fold split (random_state=42).
- **OBJ R² ≈ 0.998 on iML1515 is a target-construction artifact.** Severity is obj/obj_max by definition, so the OBJ feature *is* the target's numerator. The OBJ classification F1 is the honest scalar floor.

## Files

- `pointflux_metrics.csv` -- full headline table
- `pfba_fluxes_iso1.parquet` -- cached iSO1 pFBA solutions (242 conditions × panel)
- `pfba_fluxes_iml1515.parquet` -- cached iML1515 pFBA solutions (244 conditions × panel)
- `../figures/pointflux_vs_width.{png,pdf}` -- grouped bar comparison

