## XGBoost + SHAP 실행 커맨드 (레짐 분류 / severity 회귀)

아래 커맨드들은 **서버 기준**입니다. (프로젝트 루트에서 실행)

공통 전제:

- `results/regime_dataset.parquet`가 존재해야 합니다.
- `PYTHONPATH=$PWD/acetate_xai/src`로 소스 경로를 잡습니다.
- 출력은 **항상 새 폴더(outdir) 아래**에만 생성됩니다.

---

### 1) 레짐(label) 분류: XGBoost + SHAP

기본 실행:

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/train_xgb_shap_regime.py \
  --data results/regime_dataset.parquet \
  --outdir results/xai_xgb
```

#### 출력 폴더 / 생성 파일

- **출력 폴더**: `results/xai_xgb/`
- **생성 파일**:
  - `confusion_matrix.csv`
  - `classification_report.txt`
  - `label_mapping.csv`
  - `model.json`
  - `shap_beeswarm.png`
  - `shap_bar.png`
  - `shap_dependence_<feature>.png` (상위 3개 feature)
  - `run_metadata.json`

---

### 2) maintenance severity 회귀: XGBoost + SHAP

기본 실행:

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/train_xgb_shap_severity.py \
  --data results/regime_dataset.parquet \
  --outdir results/xai_xgb_maintenance
```

#### 출력 폴더 / 생성 파일

- **출력 폴더**: `results/xai_xgb_maintenance/`
- **생성 파일**:
  - `regression_metrics.csv` (r2, mae, rmse)
  - `model.json`
  - `shap_beeswarm.png`
  - `shap_bar.png`
  - `shap_dependence_<feature>.png` (상위 3개 feature)
  - `run_metadata.json`

---

### (선택) 빠른 스모크 테스트 옵션

서버에서 “일단 돌아가기만 하는지” 빠르게 확인하려면 다음처럼 줄여서 실행하세요.

#### 분류 스모크

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/train_xgb_shap_regime.py \
  --data results/regime_dataset.parquet \
  --outdir results/xai_xgb_smoke \
  --topk 10 \
  --n-estimators 200
```

#### 회귀 스모크

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/train_xgb_shap_severity.py \
  --data results/regime_dataset.parquet \
  --outdir results/xai_xgb_maintenance_smoke \
  --topk 10 \
  --n-estimators 200
```

