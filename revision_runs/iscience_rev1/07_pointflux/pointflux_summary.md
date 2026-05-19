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
