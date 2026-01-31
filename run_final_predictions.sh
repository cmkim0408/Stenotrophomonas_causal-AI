#!/bin/bash
# Single-command execution for server.
# Runs holdout predictions + sensitivity analysis + archive + checksums.
# NO retraining. NO model modification.

set -e
export PYTHONPATH="${PYTHONPATH:-}:$(pwd)/acetate_xai/src"

echo "[$(date -Iseconds)] Starting final predictions (holdout + sensitivity)"
python run_holdout_predictions.py
python run_sensitivity.py
echo "[$(date -Iseconds)] Predictions done"

echo "[$(date -Iseconds)] Creating archive stenotrophomonas_final_archive.tar.gz"
tar -czf stenotrophomonas_final_archive.tar.gz \
  results/holdout_predictions.csv \
  results/sensitivity \
  data/holdout_10_conditions.csv \
  2>/dev/null || true
[ -f stenotrophomonas_final_archive.tar.gz ] && echo "Archive created" || echo "Archive skipped (some paths missing)"

echo "[$(date -Iseconds)] Writing SHA256 checksums to results/CHECKSUMS.sha256"
mkdir -p results
(
  echo "# SHA256 checksums - $(date -Iseconds)"
  [ -f results/holdout_predictions.csv ] && sha256sum results/holdout_predictions.csv
  [ -f stenotrophomonas_final_archive.tar.gz ] && sha256sum stenotrophomonas_final_archive.tar.gz
  find results/sensitivity -type f -name "*.csv" -o -name "*.json" 2>/dev/null | while read f; do sha256sum "$f"; done
) > results/CHECKSUMS.sha256 2>/dev/null || true

echo "[$(date -Iseconds)] Done"
