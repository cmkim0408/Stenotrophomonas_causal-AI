## run_fva_batch_rxnfix: 여러 rxn-fix 예시 커맨드

### Bash (Linux 서버/WSL)

```bash
cd Stenotrophomonas_causal-AI

PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_fva_batch_rxnfix.py \
  --model acetate_xai/models/model.xml \
  --conditions acetate_xai/data/conditions_experiment.csv \
  --targets acetate_xai/configs/targets_120.json \
  --medium acetate_xai/configs/medium_base.yaml \
  --outdir results/test_rxnfix_multi \
  --n-jobs 1 \
  --condition-ids YE_00 \
  --rxn-fix ATPM=0 \
  --rxn-fix ATPM=20
```

### PowerShell (Windows)

```powershell
cd "C:\Cursor\Stenotrophomonas-causal AI"
$env:PYTHONPATH = "$PWD\acetate_xai\src"

python acetate_xai\scripts\run_fva_batch_rxnfix.py `
  --model acetate_xai\models\model.xml `
  --conditions acetate_xai\data\conditions_experiment.csv `
  --targets acetate_xai\configs\targets_120.json `
  --medium acetate_xai\configs\medium_base.yaml `
  --outdir results\test_rxnfix_multi `
  --n-jobs 1 `
  --condition-ids YE_00 `
  --rxn-fix ATPM=0 `
  --rxn-fix ATPM=20
```

