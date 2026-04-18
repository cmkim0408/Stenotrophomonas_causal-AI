# Feature panel ablation — summary

- N conditions: 242  ·  classes: ['N_limited', 'Ac_limited', 'O2_limited']
- Universe: 120 `width__` columns (parquet alphabetically truncated at 'FACOAL161'; M/I/N-prefixed paper anchors absent)
- Random-30 controls: 10 seeds

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
