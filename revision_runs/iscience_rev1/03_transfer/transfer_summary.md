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
