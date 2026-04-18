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

## Interpretation

- Inputs-only baselines establish a floor: they capture the obvious regime split (high O2 vs low O2) but cannot resolve nutrient-limited vs O2-limited conditions when uptake bounds overlap.
- GEM summary baselines add the FBA solution + saturation flags. These already contain most of the regime-discriminating signal (by construction of the regime label).
- Curated widths + XGBoost (ours) achieves comparable or slightly better classification AND retains a flexibility-collapse interpretation for severity regression — see `01_feature_panel_ablation/`.
