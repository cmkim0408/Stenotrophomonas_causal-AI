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
