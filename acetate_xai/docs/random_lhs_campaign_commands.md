# Random/LHS campaign (server)

This document shows recommended server commands for the large LHS campaign runner.

## Core (FBA + regime labels only; recommended default)

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_random_lhs_campaign_v2.py \
  --model acetate_xai/models/model.xml \
  --outdir results/campaigns/C_random_LHS \
  --n-samples 2000 \
  --seed 42 \
  --n-jobs 16 \
  --do-fva 0
```

Outputs:
- `results/campaigns/C_random_LHS/design.parquet`
- `results/campaigns/C_random_LHS/regime_labels.parquet`
- `results/campaigns/C_random_LHS/features.parquet` (FBA-only baseline columns)
- `results/campaigns/C_random_LHS/run_metadata.json`

## Option 1) FVA for a random subset (e.g., 500 samples)

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_random_lhs_campaign_v2.py \
  --model acetate_xai/models/model.xml \
  --outdir results/campaigns/C_random_LHS \
  --n-samples 2000 \
  --seed 42 \
  --n-jobs 16 \
  --do-fva 1 \
  --targets acetate_xai/configs/targets_120.json \
  --fraction 0.95 \
  --fva-subsample 500
```

## Option 2) FVA for all samples, but with fewer targets (e.g., 30 targets)

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_random_lhs_campaign_v2.py \
  --model acetate_xai/models/model.xml \
  --outdir results/campaigns/C_random_LHS \
  --n-samples 2000 \
  --seed 42 \
  --n-jobs 16 \
  --do-fva 1 \
  --targets acetate_xai/configs/targets_30.json \
  --fraction 0.95
```

## Option 3) FVA near phase boundaries (recommended)

This selects samples with high kNN label-disagreement scores (dense around transitions).

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_random_lhs_campaign_v2.py \
  --model acetate_xai/models/model.xml \
  --outdir results/campaigns/C_random_LHS \
  --n-samples 2000 \
  --seed 42 \
  --n-jobs 16 \
  --do-fva 1 \
  --targets acetate_xai/configs/targets_120.json \
  --fraction 0.95 \
  --fva-phase-boundary 500 \
  --fva-boundary-k 12
```

