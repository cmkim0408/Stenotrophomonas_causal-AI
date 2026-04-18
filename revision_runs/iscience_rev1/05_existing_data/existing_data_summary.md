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
