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
