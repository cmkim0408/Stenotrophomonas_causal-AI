================================================================================
FINAL PREPARATION - Holdout Predictions & Sensitivity Analysis
================================================================================
Date/Time: 2026-01-25 (UTC)

STATEMENT:
  Predictions will be generated prior to experimental validation.
  This ensures prospective validation and supports the claim:
  "Predictions were generated prior to experimental validation."

EXPLICIT NOTE:
  No retraining was performed.
  No model structure modification.
  No LHS resampling.
  No deletion of existing results.

PURPOSE:
  Code prepared for single-command execution on remote server.
  Server execution will produce timestamped logs and outputs,
  providing auditable evidence for the validation workflow.

FILES PREPARED:
  - data/holdout_10_conditions.csv    (10 holdout conditions: condition_id, acetate_mM, NH4_gL, YE_gL, pH, o2_level, o2lb)
  - run_holdout_predictions.py        (prediction-only script)
  - run_sensitivity.py                (sensitivity: frac=0.95, eps=0.5x, 2x)
  - run_final_predictions.sh          (single-command wrapper)
  - README_FINAL.txt                  (this file)

EXACT COMMANDS (run from project root on server):
  # Option A: All at once
  PYTHONPATH=$PWD/acetate_xai/src bash run_final_predictions.sh

  # Option B: Individually
  PYTHONPATH=$PWD/acetate_xai/src python run_holdout_predictions.py
  PYTHONPATH=$PWD/acetate_xai/src python run_sensitivity.py

PREREQUISITES (enforced at start - script exits if any missing):
  - results/xai_xgb/model.json (trained regime classifier)
  - results/xai_xgb_maintenance/model.json (trained severity regressor)
  - results/xai_xgb/label_mapping.csv
  (Run train_xgb_shap_regime.py and train_xgb_shap_severity.py first if missing)

HOLDOUT CSV FORMAT (data/holdout_10_conditions.csv):
  condition_id, acetate_mM, NH4_gL, YE_gL, pH, o2_level, o2lb
  - o2lb: direct model input (EX_o2_e.lb); -100=high, -50=mid, -10=low
  - Model constraints: EX_ac_e.lb=-k_ac*acetate_mM, EX_nh4_e.lb=-k_nh4*NH4_gL
  - Mapping logs printed per condition for verification

OUTPUTS (server run):
  - results/holdout_predictions.csv   (condition_id, predicted_regime, predicted_severity, git_commit, timestamp_utc, model_path)
  - results/sensitivity/frac_0.95_eps_0.5x/, frac_0.95_eps_2.0x/
  - results/sensitivity/sensitivity_summary.csv  (baseline vs eps delta)
  - results/CHECKSUMS.sha256
  - stenotrophomonas_final_archive.tar.gz

LOCAL VALIDATION (after OD available):
  cp data/holdout_od_results_template.csv data/holdout_od_results.csv
  # Fill od600_32h column
  python local_validate_holdout.py --pred results/holdout_predictions.csv --od data/holdout_od_results.csv --out local_validation
  # Outputs: local_validation/summary_metrics.csv, severity_vs_od.png, rank_plot.png, confusion_matrix.png

SERVER COMMANDS (verification sequence):
  wc -l data/holdout_10_conditions.csv          # expect 11 (header+10)
  sed -n '1,12p' data/holdout_10_conditions.csv
  ls -lh results/xai_xgb/model.json results/xai_xgb/label_mapping.csv results/xai_xgb_maintenance/model.json
  PYTHONPATH=$PWD/acetate_xai/src bash run_final_predictions.sh
  ls -lh results/holdout_predictions.csv stenotrophomonas_final_archive.tar.gz
  find results/sensitivity -maxdepth 3 -type f | head

================================================================================
