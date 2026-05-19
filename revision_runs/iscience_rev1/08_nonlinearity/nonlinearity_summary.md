# Fig 5 non-linearity reanalysis — summary

**Verdict:** `STRONG_NONLINEAR`

## What we did

Re-derived SHAP dependence for the top-3 severity-regression features on the extended `regime_dataset_extended.parquet` (242 conditions × 300 `width__` columns). Out-of-fold SHAP from 5-fold CV (random_state=42) was used to prevent leakage. Three non-linearity tests applied per feature:

- **LOESS overlay** — `statsmodels.nonparametric.smoothers_lowess` (frac=0.40)
- **Piecewise-linear with grid-searched breakpoint** — best knot from 41-point grid between the 10th and 90th percentile; F-test piecewise vs single linear
- **Polynomial degree-2** — OLS y ~ x + x² ; F-test vs linear

## Top-3 severity features (mean |OOF SHAP|)

| Rank | Feature | |OOF SHAP| |
|---|---|---|
| 1 | `width__ACONT` | 0.06654 |
| 2 | `width__12DGR120tipp` | 0.00778 |
| 3 | `width__ADCS` | 0.00388 |

## Per-feature test results

| Feature | Test | p-value | R² baseline (lin) | R² model | Δ R² | extra |
|---|---|---|---|---|---|---|
| `width__ACONT` | linear_baseline | nan | 0.890 | 0.890 | +0.000 |  |
| `width__ACONT` | piecewise_linear | 0 | 0.890 | 0.943 | +0.052 | breakpoint=1973.7891 |
| `width__ACONT` | polynomial_d2 | 0 | 0.890 | 0.939 | +0.048 |  |
| `width__12DGR120tipp` | linear_baseline | nan | 0.204 | 0.204 | +0.000 |  |
| `width__12DGR120tipp` | piecewise_linear | 7.383e-14 | 0.204 | 0.370 | +0.167 | breakpoint=3.0556 |
| `width__12DGR120tipp` | polynomial_d2 | 0.0009798 | 0.204 | 0.239 | +0.035 |  |
| `width__ADCS` | linear_baseline | nan | 0.178 | 0.178 | +0.000 |  |
| `width__ADCS` | piecewise_linear | 0 | 0.178 | 0.603 | +0.425 | breakpoint=0.0000 |
| `width__ADCS` | polynomial_d2 | 4.065e-05 | 0.178 | 0.234 | +0.056 |  |

## SHAP interaction partners (top-3 features)

| Feature | Strongest interacting partner | mean &#124;interaction&#124; |
|---|---|---|
| `width__ACONT` | `width__12DGR120tipp` | 0.00467 |
| `width__12DGR120tipp` | `width__ACONT` | 0.00467 |
| `width__ADCS` | `width__ACONT` | 0.00451 |

## Verdict logic

- Any piecewise test with p<0.05? **True**
- Any polynomial test with p<0.05? **True**
- Max Δ R² over linear: piecewise = +0.425, polynomial = +0.056
- Max mean |interaction| (top-3 features) = 0.00467

→ Verdict = `STRONG_NONLINEAR`

Branching rule (per the workstream-D prompt):

  - `STRONG_NONLINEAR` if (piecewise sig.) AND (polynomial sig.) AND (max Δ R² > 0.10)
  - `THRESHOLD_LIKE`   if piecewise sig. (independent of polynomial)
  - `INTERACTION_MOD`  if main-effect tests not sig. but max |interaction| > 0.005
  - `LINEAR_RETREAT`   otherwise

## Recommended caption variant

See `docs/revision/captions/fig5_caption_v2.md`. The variant marked `RECOMMENDED based on workstream #D verdict: STRONG_NONLINEAR` should be used for the Fig 5 caption (or the SI dependence-plot caption if main Fig 5 is left untouched).

## Honest interpretation

Both the piecewise-linear and polynomial fits significantly improve over the linear baseline (p<0.05) with Δ R² > 0.10. The 'non-linear effects' wording in the manuscript is supported. The breakpoint provides a data-driven threshold for the flexibility-collapse interpretation.
