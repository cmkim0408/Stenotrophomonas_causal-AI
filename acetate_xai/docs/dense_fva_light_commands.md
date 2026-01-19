## Dense FVA (Light) — Tier A (Meso-scale)

목표: Random LHS feasible 샘플(예: `status=="optimal"` 약 1526개) 전부에 대해 **Light FVA(30 targets)**를 수행해서 `features_dense.parquet`를 생성합니다.  
중간에 서버 세션이 끊겨도 **이미 생성된 chunk 결과는 자동 스킵**되어 이어서 실행됩니다.

---

## 입력/출력 (고정)

- **입력 모델**: `acetate_xai/models/model.xml`
- **입력 레이블**: `results/campaigns/C_random_LHS/regime_labels.parquet`
- **출력 루트**: `results/campaigns/C_random_LHS_FVA_LIGHT/`
  - `fva_parts/part_{start}_{end}.parquet`
  - `features_dense.parquet`
  - `targets_30_resolved.json`
  - `run_metadata.json`
  - `failed_samples.csv` (있을 때만)

---

## 실행 (기본)

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_dense_fva_light.py \
  --model acetate_xai/models/model.xml \
  --labels results/campaigns/C_random_LHS/regime_labels.parquet \
  --outdir results/campaigns/C_random_LHS_FVA_LIGHT \
  --chunk-size 50 \
  --n-jobs 8 \
  --fraction 0.95
```

---

## 재시작(Resume)

- 이미 생성된 `results/campaigns/C_random_LHS_FVA_LIGHT/fva_parts/part_{start}_{end}.parquet`가 있으면 해당 chunk는 자동으로 skip됩니다.
- 실패 샘플은 `failed_samples.csv`에 누적 기록되며, 전체 실행은 계속 진행됩니다.

---

## Merge만 다시 (빠르게)

chunk 결과가 이미 있고, 최종 `features_dense.parquet`만 다시 만들고 싶을 때:

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_dense_fva_light.py \
  --model acetate_xai/models/model.xml \
  --labels results/campaigns/C_random_LHS/regime_labels.parquet \
  --outdir results/campaigns/C_random_LHS_FVA_LIGHT \
  --merge-only
```

---

## 생성되는 feature 컬럼 규칙

- wide feature matrix 컬럼:
  - `width__RXNID`
  - `mid__RXNID`
- 기본 메타 컬럼:
  - `sample_id, acetate_mM, o2_uptake_max, nh4_uptake_max, atpm, objective_value, primary_regime`

