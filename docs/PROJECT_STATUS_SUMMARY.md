## 1. Why we did this (배경: BMC 리뷰 포인트 → 플랫폼 방향 전환)

본 프로젝트는 “acetate-fed 조건에서 성장(OD@32h)이 왜 비효율적인가?”라는 질문을 **GSM 기반 디지털 트윈 + 대규모 in silico intervention atlas**로 정리하고, 그 결과를 **설명 가능한 AI(XAI)**로 압축하여 논문/발표 가능한 형태의 **레짐(regime) 해석 프레임**을 만드는 것을 목표로 했습니다.  BMC/리뷰 관점에서 핵심은 “단일 in silico 결과 나열”이 아니라, **실험 endpoint(OD@32h)에 앵커링된 플랫폼형 분석**(캠페인/런/robustness/재현성)을 제시하는 것입니다.

중요하게, 본 문서에서 “causal”은 **in vivo 인과 단정**이 아니라, **in silico intervention space에서의 causal hypothesis graph(skeleton/DAG)**를 의미합니다.

## 2. Experimental anchors (OD@32h, 22조건 요약)

- **관측 데이터**: OD600 at 32 h (endpoint)
- **조건 수**: 22 conditions (1 row = 1 condition)
- **조건 축**: YE gradient / pH×YE toggle / NH4Cl gradient / acetate gradient
- **조건표 파일**: `acetate_xai/data/conditions_experiment.csv`

## 3. Digital twin & simulations (GSM, targeted FVA 120/300, campaign runs)

파이프라인은 다음 순서로 구성됩니다.

1) GSM(SBML) 로드 → 2) condition→medium(exchanges) 적용 → 3) FBA → 4) targeted FVA(주로 120 targets; 비교로 300 targets) → 5) FVA-derived feature table 생성

캠페인 산출물은 `results/campaigns/<campaign>/<run>/` 아래에 누적되며, run별로 `medium.yaml`, `fva_all.parquet`, `features.parquet`, `regime_fba.parquet`, `xai/` 등이 저장됩니다.

## 4. XAI baseline (tree/linear) 요약)

Baseline XAI는 table-first로 다음을 생성합니다.

- `results/xai/regime_table.csv`: saturation(또는 제약) 기반 레짐 요약 + narrow/wide reactions
- `results/xai/lasso_coefficients.csv`: 선형(ElasticNet/Lasso) 기반 top coefficients
- `results/xai/tree_rules.txt`: depth-limited decision tree rules

## 5. XGBoost + SHAP (classification + severity regression) 요약

### Regime classification (XGBClassifier + SHAP)

- 입력: `results/regime_dataset.parquet`의 wide feature columns(`width__/mid__/signchange__`)
- 출력: `results/xai_xgb/` (beeswarm/bar/importance + confusion matrix/report)

### Maintenance severity regression (XGBRegressor + SHAP)

- severity 정의: run_id별 max objective로 정규화한 `objective_value / max_objective_in_same_run` (0–1)
- 출력: `results/xai_xgb_maintenance/` (metrics + beeswarm/bar)

## 6. Robustness & resolution validation (0.95/0.99, 120/300)

Robustness는 다음 두 축을 중심으로 확인합니다.

- **fraction_of_optimum**: 0.95 vs 0.99 (FVA 해상도/제약감도 변화)
- **targets resolution**: 120 vs 300 targets (feature coverage/모듈 안정성)

### Robustness summary (from results/platform_summary/robustness_summary.txt)

```
Platform Robustness Summary
===========================

- Total rows: 242
- Total runs: 11
- Objective (mean/min): 0.873526 / 0.223327
- Maintenance severity (mean): 0.935883

Regime coverage (counts):
- O2_limited: 211
- N_limited: 21
- Ac_limited: 10

Interpretation note:
- 'primary_regime' is the mode of labels within each run_id.
- 'severity' is normalized within run_id by objective_value / max(objective_value).
- 'signature_modules' summarizes module frequencies derived from top linear features.
```

## 7. Causal layer (PC + bootstrap skeleton, what it means / what we claim)

본 단계는 “원인 단정”이 아니라 다음을 목표로 합니다.

- **외생 개입 변수(캠페인 토큰 파싱)** + **SHAP top-k 내부 feature** + **outcomes(primary_regime, maintenance_severity)**로 구성된 변수 집합에서
- PC Algorithm + bootstrap으로 **edge stability(빈도)**를 계산하여
- **in silico causal hypothesis graph**를 제시합니다.

이 그래프는 “검증 가능한 가설”을 제공하며, in vivo causal claim을 대체하지 않습니다.

## 8. Where to find results (핵심 파일 경로 목록)

아래는 이 repo의 결과 패키지(`results_all_v2.tar.gz`)를 풀었을 때의 `results/` 기준 핵심 경로입니다.

- `results/platform_summary/`
  - `platform_summary.csv`, `regime_map.png`, `robustness_summary.txt`, `signature_modules.csv`, `regime_coverage.csv`
- `results/xai_xgb/`
  - `shap_beeswarm.png`, `shap_bar.png`, `shap_importance.csv`, `confusion_matrix.csv`, `classification_report.txt`
- `results/xai_xgb_maintenance/`
  - `regression_metrics.csv`, `shap_beeswarm.png`, `shap_bar.png`
- `results/causal_dag/`
  - `causal_dag.png`/`causal_dag.pdf`, `edge_stability.csv`, `dag_edges.csv`
- `results/regime_dataset.parquet`
- 캠페인 구조:
  - `results/campaigns/<campaign>/<run>/{medium.yaml, features.parquet, fva_all.parquet, regime_fba.parquet, xai/}`

## 9. Suggested paper figure list (5장) + 어떤 파일이 근거인지

1) **Workflow overview**: `results/figures_draft/Fig01_workflow.png`
2) **Regime map**: `results/platform_summary/regime_map.png` → `results/figures_draft/Fig02_regime_map.png`
3) **Regime classifier SHAP**: `results/xai_xgb/shap_beeswarm.png` → `results/figures_draft/Fig03_shap_regime_beeswarm.png`
4) **Severity regressor SHAP**: `results/xai_xgb_maintenance/shap_beeswarm.png` → `results/figures_draft/Fig04_shap_severity_beeswarm.png`
5) **Causal hypothesis skeleton**: `results/causal_dag/causal_dag.png` → `results/figures_draft/Fig05_causal_dag.png` + `edge_stability.csv` 기반 `Fig06_top_edges_bar.png`

