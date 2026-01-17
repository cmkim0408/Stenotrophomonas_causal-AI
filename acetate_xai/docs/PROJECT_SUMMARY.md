## 1. Background and Motivation

이 프로젝트의 목표는 **Stenotrophomonas SO-1의 acetate 기반 성장 비효율(낮은 OD@32h 패턴)**을 **Genome-scale metabolic model(GSM) 기반 디지털 트윈**으로 재현하고, 대량 in-silico 조건에서 생성된 대사적 제약/유연성 신호를 **설명 가능한 AI(XAI)**로 요약해 **조건 의존적 성장 제한 레짐(regime)**을 정리하는 것입니다.  
중요하게, 본 프로젝트의 “causal”은 **in vivo 인과 단정**이 아니라, **in silico intervention space에서의 인과 가설 그래프(hypothesis DAG)** 생성/검증을 의미합니다.

### Abstract-friendly one-paragraph summary (no overclaim)

We constructed a genome-scale metabolic digital twin of *Stenotrophomonas* SO-1 and generated a large intervention-driven simulation atlas across controlled media perturbations. Using targeted flux variability signatures and mechanistic constraint indicators, we built explainable predictors of OD600 endpoints (32 h) and mapped condition-dependent growth-limiting regimes. Finally, we derived an in silico causal hypothesis graph (PC + bootstrap stability) that links external interventions to internal metabolic flexibility features and outcomes, providing testable mechanistic hypotheses without claiming in vivo causality.

---

## 2. Experimental Dataset (OD at 32 h)

### What we have

- **Measured endpoint**: **OD600 at 32 h**  
- **Total conditions**: **22** (long-format, 1 row = 1 condition)
- **Stored at**: `acetate_xai/data/conditions_experiment.csv`

### Condition sets (conceptual)

- **Yeast extract gradient**: yeast_extract_gL sweep (pH0=7.0, others fixed)
- **pH × YE toggle**: pH0 in {6.0, 6.5, 7.0, 7.5} × YE present/absent
- **NH4Cl gradient**: nh4cl_gL sweep (YE fixed to 0.2 g/L, pH0=7.0)
- **Acetate gradient**: acetate_mM sweep (pH0=7.0, others fixed)

---

## 3. Digital-Twin and In Silico Pipeline

### Core flow (what the code does)

1) **Load GSM (SBML)**  
2) **Apply base medium + condition → exchange bounds**  
3) **FBA** to obtain objective value (growth proxy)  
4) **Targeted FVA** on **120 reactions** with `fraction_of_optimum` (e.g., 0.95)  
5) **Feature engineering** from FVA min/max  
6) **XAI reports** (baseline linear/tree; plus XGBoost/SHAP)  
7) **Robustness checks** (fraction, targets, etc.)  
8) **Causal layer**: interventions + SHAP-selected features + outcomes → PC + bootstrap stability

### Key implementation scripts (high level)

- **Batch targeted FVA**: `acetate_xai/scripts/run_fva_batch.py` (+ variants)
- **Collect + build features**: `acetate_xai/scripts/collect_fva_parts.py`
- **Baseline XAI report**: `acetate_xai/scripts/xai_report.py`
- **Mechanistic regime by FBA flux+bound**: `acetate_xai/scripts/run_fba_regime.py` (+ rxnfix variant)
- **Platform dataset + classifiers**: `acetate_xai/scripts/build_regime_dataset.py`, `acetate_xai/scripts/train_regime_xai.py`
- **Maintenance severity regression**: `acetate_xai/scripts/train_maintenance_xai.py`
- **XGB+SHAP**: `acetate_xai/scripts/train_xgb_shap_regime.py`, `acetate_xai/scripts/train_xgb_shap_severity.py`
- **Causal**: `acetate_xai/scripts/build_causal_dataset.py`, `acetate_xai/scripts/export_shap_topk.py`, `acetate_xai/scripts/run_causal_discovery.py`, `acetate_xai/scripts/plot_causal_dag.py`

---

## 4. Campaigns and Runs (What was varied)

캠페인은 `results/campaigns/<campaign>/<run>/...` 구조로 저장됩니다. 각 run은 보통 다음을 포함합니다:

- `medium.yaml`: 해당 run에서 사용한 medium/config 스냅샷
- `fva_parts/condition_id=*.parquet`: 조건별 분산 FVA 결과(120 reactions)
- `fva_all.parquet`: 집계(long)
- `features.parquet`: wide feature matrix (XAI 입력)
- `regime_fba.parquet`: FBA 기반 포화/제약 레짐 테이블
- `xai/`: baseline XAI outputs

Release(`results_all_v2.tar.gz`)에서 확인된 캠페인 예시는 다음과 같습니다.

- `C1_o2_sweep`: oxygen lower bound sweep
- `C2_nh4_free`: NH4 cap loosen / mode change
- `C3_o2_mid`: intermediate oxygen caps
- `C3_o2_mid_300targets`: targets 수 변화(예: 300 targets)
- `C3_o2_mid_f099`: fraction_of_optimum 변화(예: 0.99)
- `C4_atpm_sweep`, `C4_atpm_sweep_v2`: ATPM 고정/스윕 (v2는 rxnfix 기반 산출물 우선)
- `C5_pi_sweep`: phosphate lower bound sweep

---

## 5. Explainable AI (Baseline + XGBoost/SHAP)

### Baseline XAI (table-first)

- 입력: `results/features.parquet`
- 출력: `results/xai/` 아래
  - `regime_table.csv`
  - `lasso_coefficients.csv` (ElasticNet 기반 top coefficients)
  - `tree_rules.txt` (depth-limited tree)

### XGBoost/SHAP (stronger nonlinear model + global/local explanation)

분류(레짐):

- `acetate_xai/scripts/train_xgb_shap_regime.py`
- 출력: `results/xai_xgb/`
  - `shap_beeswarm.png`, `shap_bar.png`, `shap_importance.csv`
  - `confusion_matrix.csv`, `classification_report.txt`, `label_mapping.csv`
  - `model.json`, `run_metadata.json`, `shap_dependence_*.png`

회귀(maintenance severity):

- `acetate_xai/scripts/train_xgb_shap_severity.py`
- 출력: `results/xai_xgb_maintenance/`
  - `regression_metrics.csv`
  - `shap_beeswarm.png`, `shap_bar.png`
  - `model.json`, `run_metadata.json`, `shap_dependence_*.png`

관련 실행 커맨드 문서:

- `acetate_xai/docs/xgb_shap_commands.md`

---

## 6. Robustness and Resolution Validation (0.95 vs 0.99, 120 vs 300)

본 프로젝트의 “robustness”는 크게 두 축을 통해 확보합니다.

- **Optimization fraction**: `fraction_of_optimum` (예: 0.95 vs 0.99)  
  - 준최적 허용 시 FVA 폭이 달라지고, OD endpoint의 변동성을 더 잘 반영할 수 있음
- **Target resolution**: targeted FVA reaction 수 (예: 120 vs 300 targets)  
  - feature 해상도/모듈 커버리지를 바꾸며, XAI top features가 안정적인지 확인

플랫폼 요약 산출물에서 전체적으로 정리됩니다:

- `results/platform_summary/robustness_summary.txt`
- `results/platform_summary/signature_modules.csv`

---

## 7. Causal Layer (PC + bootstrap; skeleton graph)

목표는 “원인 단정”이 아니라 **in silico causal hypothesis graph**를 생성하는 것입니다.

### Variable construction (high level)

- **Exogenous interventions** (가능하면 항상 포함):
  - `atpm_fixed`, `o2_lb`, `acetate_mode`, `nh4_mode`, `frac_opt`, `targets_n`
- **SHAP top-k features**:
  - regime top-k + severity top-k (중복 제거, 필요 시 자동 축소)
- **Outcomes**:
  - `primary_regime`
  - `maintenance_severity`

### Discovery + stability

- 방법: PC Algorithm (`fisherz` 기본; 필요 시 `kci` 가능)  
- Bootstrap: edge stability(빈도) 계산  
- 출력: edge list + stability table + metadata JSON  
- 시각화: networkx + matplotlib

---

## 8. Key Output Files and Where to Find Them

아래 경로는 **Release 압축(`results_all_v2.tar.gz`)을 풀었을 때의 `results/` 기준**입니다.

### Platform summary (final “what happened”)

- `results/platform_summary/`
  - `platform_summary.csv`: run_id별 primary_regime, objective 통계, severity 통계 요약
  - `regime_map.png`: primary_regime vs maintenance_severity scatter (run_id 색상)
  - `robustness_summary.txt`: 수치 삽입된 자동 요약 텍스트
  - `signature_modules.csv`: 분류/회귀 top feature에서 모듈 빈도 집계
  - `regime_coverage.csv`: 전체 레짐 분포(커버리지)

### XGBoost + SHAP (regime classifier)

- `results/xai_xgb/`
  - `shap_beeswarm.png`: SHAP summary (beeswarm)
  - `shap_bar.png`: SHAP global bar
  - `shap_importance.csv`: (feature, mean_abs_shap) 정렬 테이블
  - `confusion_matrix.csv`: 분류 혼동행렬
  - `classification_report.txt`: precision/recall/f1

### XGBoost + SHAP (maintenance severity regressor)

- `results/xai_xgb_maintenance/`
  - `regression_metrics.csv`: r2/mae/rmse
  - `shap_beeswarm.png`: SHAP summary (beeswarm)
  - `shap_bar.png`: SHAP global bar

### Causal DAG outputs

- `results/causal_dag/`
  - `causal_dag.png` / `causal_dag.pdf`: DAG visualization (edge thickness/alpha = stability)
  - `edge_stability.csv`: bootstrap edge frequency
  - `dag_edges.csv`: discovered edges

### Core merged datasets

- `results/regime_dataset.parquet`: 캠페인 전반 run들을 통합한 레짐 플랫폼 데이터셋

### Campaign folder structure (per run)

- `results/campaigns/<campaign>/<run>/{medium.yaml, features.parquet, fva_all.parquet, regime_fba.parquet, xai/}`

---

## 9. Reproducibility: How to Re-run (Server commands)

아래는 **서버에서 재현 가능한 “정석 순서”**입니다. (프로젝트 루트에서 실행)

### A) 설치/환경

```bash
git clone https://github.com/cmkim0408/Stenotrophomonas_causal-AI.git
cd Stenotrophomonas_causal-AI/acetate_xai
python3 -m venv .venv
source .venv/bin/activate
pip install -e ".[dev]"
```

### B) (이미 결과가 있다면) 플랫폼 요약/지도 생성

```bash
PYTHONPATH=$PWD/src python scripts/summarize_platform_results.py \
  --dataset results/regime_dataset.parquet \
  --xai-dir results/xai_platform \
  --maintenance-dir results/xai_platform_maintenance \
  --outdir results/platform_summary

PYTHONPATH=$PWD/src python scripts/plot_regime_map.py \
  --dataset results/regime_dataset.parquet \
  --outdir results/platform_summary
```

### C) XGBoost/SHAP

```bash
PYTHONPATH=$PWD/src python scripts/train_xgb_shap_regime.py \
  --data results/regime_dataset.parquet \
  --outdir results/xai_xgb

PYTHONPATH=$PWD/src python scripts/train_xgb_shap_severity.py \
  --data results/regime_dataset.parquet \
  --outdir results/xai_xgb_maintenance
```

### D) Causal layer (dataset → topk → PC + bootstrap → plot)

```bash
PYTHONPATH=$PWD/src python scripts/build_causal_dataset.py \
  --in results/regime_dataset.parquet \
  --out results/causal_dataset.parquet

PYTHONPATH=$PWD/src python scripts/export_shap_topk.py \
  --data results/causal_dataset.parquet \
  --outdir results/causal_topk \
  --topk 15

PYTHONPATH=$PWD/src python scripts/run_causal_discovery.py \
  --data results/causal_dataset.parquet \
  --topk-regime results/causal_topk/top_features_regime.csv \
  --topk-severity results/causal_topk/top_features_severity.csv \
  --outdir results/causal_dag \
  --bootstrap 50 \
  --method pc

PYTHONPATH=$PWD/src python scripts/plot_causal_dag.py \
  --edges results/causal_dag/dag_edges.csv \
  --stability results/causal_dag/edge_stability.csv \
  --outdir results/causal_dag
```

---

## 10. Public Release Artifacts (GitHub Release)

프로젝트 배포는 2축으로 구성됩니다.

- **GitHub Repo**: 코드 + 설정 + (필요시) 모델 파일  
- **GitHub Release**: 대량 결과 산출물 패키지

Release 정보:

- **tag**: `results-2026-01-16`
- **latest archive**: `results_all_v2.tar.gz`
- **sha256**: `eb76e5cedaaad6cb575c012e5bd2b5c895c501548ec50ac51e474dafde9c89cf`

무결성 검증 예시:

```bash
sha256sum results_all_v2.tar.gz
```

---

## 11. Suggested Figures (5-6 figures) and What Each Shows

1) **Workflow schematic**: OD@32h → constraint/targeted FVA → features → XAI → causal hypothesis DAG  
2) **Regime coverage + map**: `platform_summary/regime_coverage.csv` + `regime_map.png`  
3) **Mechanistic regime table**: FBA-based saturation(acetate/o2/nh4/pi) + objective trend  
4) **XGB+SHAP (classification)**: beeswarm + bar + top features table  
5) **XGB+SHAP (severity regression)**: beeswarm + bar + metrics table  
6) **Causal hypothesis DAG**: `causal_dag.png` + `edge_stability.csv` 요약(상위 안정 edge)

---

## 12. Current Status Checklist (Done/To-do)

### Done

- [x] OD@32h 조건표 표준화(22조건) 및 파이프라인 연결
- [x] GSM 로드 → condition→medium 적용 → FBA → targeted FVA(120) 배치 실행
- [x] 분산 parquet 집계 → features.parquet 생성
- [x] baseline XAI 리포트 생성 (table-first)
- [x] 플랫폼 데이터셋(`regime_dataset.parquet`) 및 분류/회귀 요약 생성
- [x] XGBoost/SHAP 분류/회귀 결과 생성 및 파일화
- [x] in silico causal hypothesis DAG (PC + bootstrap) 및 시각화
- [x] Release 패키징(`results_all_v2.tar.gz`) + sha256 공개

### To-do (optional / next iteration)

- [ ] causal 변수/priors 정교화(도메인 지식 기반 금지/허용 edge 제약)
- [ ] independence test(혼합형 데이터) 정교화(continuous vs categorical)
- [ ] figure polishing (paper-ready aesthetics, captions)

