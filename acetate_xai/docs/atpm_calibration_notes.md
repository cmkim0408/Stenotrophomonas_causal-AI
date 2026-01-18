# ATPM Calibration Notes (Anchor-based Parameter Estimation)

## Why ATPM calibration?

Acetate는 weak acid stress를 유발할 수 있고, 이 경우 세포는 **proton/ion homeostasis**를 유지하기 위해 추가적인 **ATP maintenance demand(ATPM_eff)**를 필요로 할 수 있습니다.  
GSM/FBA에서 이 효과는 보통 `ATPM` 반응(maintenance) 플럭스를 증가시키는 형태로 모델링합니다.

본 프로젝트에서는 **실험 endpoint(OD600 at 32 h)**에 앵커링하여,
소수의 **anchor 조건**에서 `ATPM_eff`를 역추정하고,
이를 \(ATPM\_eff = f(\text{acetate\_mM})\) 형태(1차 버전: 선형)로 전체 조건에 적용한 뒤,
**objective_value vs measured_OD의 Spearman 상관(ρ)**가 개선되는지 확인합니다.

중요: 이는 **parameter estimation(캘리브레이션)**이며, in vivo 인과를 주장하지 않습니다.

## Train/Test split (anchors vs held-out)

- **Train (anchors)**: acetate_gradient 세트에서 YE=0, pH0=7.0 조건을 기본 앵커로 사용  
  기본값: `AC_25`, `AC_100`, `AC_150` (최소 2개면 linear fit 가능)
- **Test (held-out)**: 나머지 조건들(그리고 필요하면 앵커도 포함한 전체 22조건)에서 재시뮬레이션 결과의 일반화 성능을 확인

현재 구현은:
- calibration은 **anchor id만 사용**
- 재시뮬레이션(run)은 **전체 조건에 f(Ac)를 적용**
- Figure 2 Panel E before/after는 **전체 조건 기준**으로 ρ를 계산 (요구사항)

## Calibration objective (v1)

v1은 구현 단순화를 위해 각 anchor에서 다음을 최소화하는 방식으로 `atpm_best`를 선택합니다.

- \(OD\_{norm} = (OD - OD_{min})/(OD_{max}-OD_{min})\)  (anchor들 사이에서 정규화)
- \(\\mu\_{norm}(ATPM) = objective(ATPM) / objective(ATPM=0)\)
- 각 anchor에서 \(|\\mu\_{norm} - OD\_{norm}|\)가 최소가 되는 ATPM을 `atpm_best`로 선택

옵션으로 `--mode rank`를 주면 (anchor가 작을 때) rank matching 기반으로도 선택 가능합니다.

## Model application: f(Ac) → per-condition rxn-fix

캘리브레이션 결과는 `atpm_fit.json`에 저장되며, v1은 선형 모델입니다.

- \(ATPM\_{eff} = clip(a + b \\cdot acetate\\_mM, 0, 200)\)
- 각 조건에서 medium 적용 후, `ATPM` 반응을 **fixed flux(lb=ub=ATPM_eff)**로 고정한 다음 FBA/FVA를 실행합니다.

## Wording (no over-claim)

권장 표현(논문/발표용 템플릿):

- “We performed **anchor-based parameter estimation** of ATP maintenance demand to better align digital-twin growth proxies with OD endpoints.”
- “The calibrated parameterization improves the **in silico growth–OD rank agreement** under acetate stress.”
- “This calibration supports an **in silico hypothesis** that maintenance demand increases with acetate load.”

피해야 할 과장 표현:

- “We proved acetate causes maintenance increase in vivo.” (금지)
- “Causal effect is confirmed experimentally.” (금지)

## Server commands (end-to-end)

아래는 서버에서 “캘리브레이션 → 재시뮬레이션 → Figure 2E before/after 생성”을 한 번에 수행하는 최소 커맨드입니다.

### 1) Calibration (anchors only)

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/calibrate_atpm_from_anchors.py \
  --model acetate_xai/models/model.xml \
  --conditions acetate_xai/data/conditions_experiment.csv \
  --medium acetate_xai/configs/medium_base.yaml \
  --anchor-ids AC_25 AC_100 AC_150 \
  --atpm-min 0 --atpm-max 200 --atpm-step 5 \
  --mode norm \
  --out results/campaigns/C6_atpm_calibrated/calibration/
```

출력(필수):
- `results/campaigns/C6_atpm_calibrated/calibration/anchor_scan.parquet`
- `results/campaigns/C6_atpm_calibrated/calibration/anchor_best.csv`
- `results/campaigns/C6_atpm_calibrated/calibration/atpm_fit.json`

### 2) Calibrated run (all 22 conditions)

```bash
PYTHONPATH=$PWD/acetate_xai/src python acetate_xai/scripts/run_campaign_atpm_calibrated.py \
  --model acetate_xai/models/model.xml \
  --conditions acetate_xai/data/conditions_experiment.csv \
  --targets acetate_xai/configs/targets_120.json \
  --medium acetate_xai/configs/medium_base.yaml \
  --regime-config acetate_xai/configs/regime_exchanges.yaml \
  --atpm-fit results/campaigns/C6_atpm_calibrated/calibration/atpm_fit.json \
  --outdir results/campaigns/C6_atpm_calibrated/run__atpmCalib_linear/ \
  --n-jobs 8 \
  --fraction 0.95
```

출력(필수):
- `.../features.parquet`
- `.../regime_fba.parquet`
- `.../xai/regime_table.csv`

### 3) Figures (Fig02E before/after + rho_compare)

```bash
python acetate_xai/scripts/make_paper_figures.py
```

출력(필수):
- `results/figures_draft/Fig02E_before_objective_vs_OD.png`
- `results/figures_draft/Fig02E_after_calibObjective_vs_OD.png`
- `results/figures_draft/Fig02E_rho_compare.txt`

