from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np
import pandas as pd


def _find_results_root() -> Path:
    """
    Prefer extracted release results over any partial local results/.
    """
    candidates = [
        Path("results_all_v2_extracted") / "results",
        Path("results_all_extracted") / "results",
        Path("results"),
    ]
    for p in candidates:
        if p.exists() and p.is_dir():
            return p
    raise FileNotFoundError(f"No results/ directory found (tried: {', '.join(str(p) for p in candidates)})")


def _copy(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dst)


def _write_status_summary(results_root: Path, out_md: Path) -> None:
    """
    Create a 1–2 page internal summary (markdown) based on known outputs.
    """
    # Light-touch parsing (optional): robustness_summary.txt text inclusion
    robustness_path = results_root / "platform_summary" / "robustness_summary.txt"
    robustness_txt = ""
    if robustness_path.exists():
        robustness_txt = robustness_path.read_text(encoding="utf-8", errors="replace").strip()

    out_md.parent.mkdir(parents=True, exist_ok=True)

    parts: list[str] = []
    parts.append("## 1. Why we did this (배경: BMC 리뷰 포인트 → 플랫폼 방향 전환)")
    parts.append("")
    parts.append(
        "본 프로젝트는 “acetate-fed 조건에서 성장(OD@32h)이 왜 비효율적인가?”라는 질문을 **GSM 기반 디지털 트윈 + 대규모 in silico intervention atlas**로 정리하고, "
        "그 결과를 **설명 가능한 AI(XAI)**로 압축하여 논문/발표 가능한 형태의 **레짐(regime) 해석 프레임**을 만드는 것을 목표로 했습니다.  "
        "BMC/리뷰 관점에서 핵심은 “단일 in silico 결과 나열”이 아니라, **실험 endpoint(OD@32h)에 앵커링된 플랫폼형 분석**(캠페인/런/robustness/재현성)을 제시하는 것입니다."
    )
    parts.append("")
    parts.append(
        "중요하게, 본 문서에서 “causal”은 **in vivo 인과 단정**이 아니라, **in silico intervention space에서의 causal hypothesis graph(skeleton/DAG)**를 의미합니다."
    )
    parts.append("")
    parts.append("## 2. Experimental anchors (OD@32h, 22조건 요약)")
    parts.append("")
    parts.append("- **관측 데이터**: OD600 at 32 h (endpoint)")
    parts.append("- **조건 수**: 22 conditions (1 row = 1 condition)")
    parts.append("- **조건 축**: YE gradient / pH×YE toggle / NH4Cl gradient / acetate gradient")
    parts.append("- **조건표 파일**: `acetate_xai/data/conditions_experiment.csv`")
    parts.append("")
    parts.append("## 3. Digital twin & simulations (GSM, targeted FVA 120/300, campaign runs)")
    parts.append("")
    parts.append("파이프라인은 다음 순서로 구성됩니다.")
    parts.append("")
    parts.append(
        "1) GSM(SBML) 로드 → 2) condition→medium(exchanges) 적용 → 3) FBA → 4) targeted FVA(주로 120 targets; 비교로 300 targets) → 5) FVA-derived feature table 생성"
    )
    parts.append("")
    parts.append(
        "캠페인 산출물은 `results/campaigns/<campaign>/<run>/` 아래에 누적되며, run별로 `medium.yaml`, `fva_all.parquet`, `features.parquet`, `regime_fba.parquet`, `xai/` 등이 저장됩니다."
    )
    parts.append("")
    parts.append("## 4. XAI baseline (tree/linear) 요약)")
    parts.append("")
    parts.append("Baseline XAI는 table-first로 다음을 생성합니다.")
    parts.append("")
    parts.append("- `results/xai/regime_table.csv`: saturation(또는 제약) 기반 레짐 요약 + narrow/wide reactions")
    parts.append("- `results/xai/lasso_coefficients.csv`: 선형(ElasticNet/Lasso) 기반 top coefficients")
    parts.append("- `results/xai/tree_rules.txt`: depth-limited decision tree rules")
    parts.append("")
    parts.append("## 5. XGBoost + SHAP (classification + severity regression) 요약")
    parts.append("")
    parts.append("### Regime classification (XGBClassifier + SHAP)")
    parts.append("")
    parts.append("- 입력: `results/regime_dataset.parquet`의 wide feature columns(`width__/mid__/signchange__`)")
    parts.append("- 출력: `results/xai_xgb/` (beeswarm/bar/importance + confusion matrix/report)")
    parts.append("")
    parts.append("### Maintenance severity regression (XGBRegressor + SHAP)")
    parts.append("")
    parts.append("- severity 정의: run_id별 max objective로 정규화한 `objective_value / max_objective_in_same_run` (0–1)")
    parts.append("- 출력: `results/xai_xgb_maintenance/` (metrics + beeswarm/bar)")
    parts.append("")
    parts.append("## 6. Robustness & resolution validation (0.95/0.99, 120/300)")
    parts.append("")
    parts.append("Robustness는 다음 두 축을 중심으로 확인합니다.")
    parts.append("")
    parts.append("- **fraction_of_optimum**: 0.95 vs 0.99 (FVA 해상도/제약감도 변화)")
    parts.append("- **targets resolution**: 120 vs 300 targets (feature coverage/모듈 안정성)")
    parts.append("")
    if robustness_txt:
        parts.append("### Robustness summary (from results/platform_summary/robustness_summary.txt)")
        parts.append("")
        parts.append("```")
        parts.append(robustness_txt)
        parts.append("```")
        parts.append("")
    parts.append("## 7. Causal layer (PC + bootstrap skeleton, what it means / what we claim)")
    parts.append("")
    parts.append("본 단계는 “원인 단정”이 아니라 다음을 목표로 합니다.")
    parts.append("")
    parts.append("- **외생 개입 변수(캠페인 토큰 파싱)** + **SHAP top-k 내부 feature** + **outcomes(primary_regime, maintenance_severity)**로 구성된 변수 집합에서")
    parts.append("- PC Algorithm + bootstrap으로 **edge stability(빈도)**를 계산하여")
    parts.append("- **in silico causal hypothesis graph**를 제시합니다.")
    parts.append("")
    parts.append("이 그래프는 “검증 가능한 가설”을 제공하며, in vivo causal claim을 대체하지 않습니다.")
    parts.append("")
    parts.append("## 8. Where to find results (핵심 파일 경로 목록)")
    parts.append("")
    parts.append("아래는 이 repo의 결과 패키지(`results_all_v2.tar.gz`)를 풀었을 때의 `results/` 기준 핵심 경로입니다.")
    parts.append("")
    parts.append("- `results/platform_summary/`")
    parts.append("  - `platform_summary.csv`, `regime_map.png`, `robustness_summary.txt`, `signature_modules.csv`, `regime_coverage.csv`")
    parts.append("- `results/xai_xgb/`")
    parts.append("  - `shap_beeswarm.png`, `shap_bar.png`, `shap_importance.csv`, `confusion_matrix.csv`, `classification_report.txt`")
    parts.append("- `results/xai_xgb_maintenance/`")
    parts.append("  - `regression_metrics.csv`, `shap_beeswarm.png`, `shap_bar.png`")
    parts.append("- `results/causal_dag/`")
    parts.append("  - `causal_dag.png`/`causal_dag.pdf`, `edge_stability.csv`, `dag_edges.csv`")
    parts.append("- `results/regime_dataset.parquet`")
    parts.append("- 캠페인 구조:")
    parts.append("  - `results/campaigns/<campaign>/<run>/{medium.yaml, features.parquet, fva_all.parquet, regime_fba.parquet, xai/}`")
    parts.append("")
    parts.append("## 9. Suggested paper figure list (5장) + 어떤 파일이 근거인지")
    parts.append("")
    parts.append("1) **Workflow overview**: `results/figures_draft/Fig01_workflow.png`")
    parts.append("2) **Regime map**: `results/platform_summary/regime_map.png` → `results/figures_draft/Fig02_regime_map.png`")
    parts.append("3) **Regime classifier SHAP**: `results/xai_xgb/shap_beeswarm.png` → `results/figures_draft/Fig03_shap_regime_beeswarm.png`")
    parts.append("4) **Severity regressor SHAP**: `results/xai_xgb_maintenance/shap_beeswarm.png` → `results/figures_draft/Fig04_shap_severity_beeswarm.png`")
    parts.append("5) **Causal hypothesis skeleton**: `results/causal_dag/causal_dag.png` → `results/figures_draft/Fig05_causal_dag.png` + `edge_stability.csv` 기반 `Fig06_top_edges_bar.png`")
    parts.append("")

    content = "\n".join(parts) + "\n"

    out_md.write_text(content, encoding="utf-8")


def _plot_workflow(out_png: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import FancyBboxPatch

    out_png.parent.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(13, 4.2))
    ax.set_axis_off()

    boxes = [
        "Experiment anchors\n(OD@32h, 22 conditions)",
        "GSM digital twin\n(SBML model)",
        "FBA + Targeted FVA\n(120/300 targets)",
        "Feature table\n(width/mid/signchange)",
        "Regime labeling\n(saturation + rules)",
        "XAI\n(baseline + XGB/SHAP)",
        "Robustness / resolution\n(0.95 vs 0.99; 120 vs 300)",
        "Causal skeleton\n(PC + bootstrap stability)",
    ]

    # positions in axes fraction coordinates
    xs = np.linspace(0.03, 0.97, len(boxes))
    y = 0.55
    w = 0.115
    h = 0.32

    for i, (x, text) in enumerate(zip(xs, boxes)):
        x0 = x - w / 2
        rect = FancyBboxPatch(
            (x0, y),
            w,
            h,
            boxstyle="round,pad=0.02,rounding_size=0.02",
            linewidth=1.2,
            facecolor="white",
            edgecolor="black",
            transform=ax.transAxes,
        )
        ax.add_patch(rect)
        ax.text(x, y + h / 2, text, ha="center", va="center", fontsize=9, transform=ax.transAxes)

        # arrows
        if i < len(boxes) - 1:
            ax.annotate(
                "",
                xy=(x + w / 2, y + h / 2),
                xytext=(xs[i + 1] - w / 2, y + h / 2),
                arrowprops=dict(arrowstyle="->", lw=1.2),
                xycoords=ax.transAxes,
                textcoords=ax.transAxes,
            )

    ax.set_title("Acetate-fed Growth Diagnosis Platform (Experiment → Digital Twin → XAI → Causal Hypothesis)", fontsize=12)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_top_edges_bar(edge_stability_csv: Path, out_png: Path, topn: int = 10) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_png.parent.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(edge_stability_csv)
    if "edge" not in df.columns or "frequency" not in df.columns:
        raise ValueError("edge_stability.csv must contain columns: edge, frequency")
    df["frequency"] = pd.to_numeric(df["frequency"], errors="coerce").fillna(0.0)
    top = df.sort_values("frequency", ascending=False).head(int(topn)).copy()
    top = top.sort_values("frequency", ascending=True)  # horizontal bar: bottom = largest

    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.barh(top["edge"].astype(str), top["frequency"])
    ax.set_xlabel("edge stability frequency (bootstrap)")
    ax.set_ylabel("edge")
    ax.set_title("Top stable edges (PC + bootstrap)")
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> int:
    results_root = _find_results_root()

    # New outputs only
    figures_out = Path("results") / "figures_draft"
    doc_out = Path("docs") / "PROJECT_STATUS_SUMMARY.md"

    # Source files inside results_root
    src_regime_map = results_root / "platform_summary" / "regime_map.png"
    src_platform_csv = results_root / "platform_summary" / "platform_summary.csv"
    src_robustness_txt = results_root / "platform_summary" / "robustness_summary.txt"
    src_shap_regime = results_root / "xai_xgb" / "shap_beeswarm.png"
    src_shap_sev = results_root / "xai_xgb_maintenance" / "shap_beeswarm.png"
    src_causal_png = results_root / "causal_dag" / "causal_dag.png"
    src_edge_stab = results_root / "causal_dag" / "edge_stability.csv"

    required = [
        src_regime_map,
        src_platform_csv,
        src_robustness_txt,
        src_shap_regime,
        src_shap_sev,
        src_causal_png,
        src_edge_stab,
    ]
    missing = [p for p in required if not p.exists()]
    if missing:
        msg = "\n".join(f"- {p}" for p in missing)
        raise FileNotFoundError(f"Missing required input files under {results_root}:\n{msg}")

    # A) document
    _write_status_summary(results_root, doc_out)

    # B) figures
    figures_out.mkdir(parents=True, exist_ok=True)

    # Fig01: workflow (drawn)
    _plot_workflow(figures_out / "Fig01_workflow.png")

    # Fig02-Fig05: copies
    _copy(src_regime_map, figures_out / "Fig02_regime_map.png")
    _copy(src_shap_regime, figures_out / "Fig03_shap_regime_beeswarm.png")
    _copy(src_shap_sev, figures_out / "Fig04_shap_severity_beeswarm.png")
    _copy(src_causal_png, figures_out / "Fig05_causal_dag.png")

    # Fig06: top edges bar (drawn)
    _plot_top_edges_bar(src_edge_stab, figures_out / "Fig06_top_edges_bar.png", topn=10)

    print(f"[OK] Wrote document: {doc_out}")
    print(f"[OK] Wrote figures under: {figures_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

