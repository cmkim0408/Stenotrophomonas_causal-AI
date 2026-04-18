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

