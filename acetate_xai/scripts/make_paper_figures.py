from __future__ import annotations

import argparse
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

    # wider canvas + smaller boxes to increase spacing
    fig, ax = plt.subplots(figsize=(18, 4.2))
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
    xs = np.linspace(0.04, 0.96, len(boxes))
    y = 0.55
    w = 0.095
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
                # matplotlib arrow goes from xytext -> xy; force rightward arrows
                xy=(xs[i + 1] - w / 2, y + h / 2),
                xytext=(x + w / 2, y + h / 2),
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


def _plot_regime_map_primary_only(regime_map_csv: Path, out_png: Path) -> None:
    """
    Rebuild regime map without run_id legend:
      - x: primary_regime (categorical)
      - y: maintenance_severity
      - color: primary_regime
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_png.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(regime_map_csv)
    need = {"primary_regime", "maintenance_severity"}
    if not need.issubset(df.columns):
        raise ValueError(f"regime_map.csv missing required columns: {sorted(need)}")

    df["primary_regime"] = df["primary_regime"].astype(str).str.strip()
    df["maintenance_severity"] = pd.to_numeric(df["maintenance_severity"], errors="coerce")
    df = df.dropna(subset=["primary_regime", "maintenance_severity"]).copy()
    df = df[df["primary_regime"].str.lower() != "nan"]

    regimes = sorted(df["primary_regime"].unique().tolist())
    x_map = {r: i for i, r in enumerate(regimes)}
    x = df["primary_regime"].map(x_map).astype(float).values
    y = df["maintenance_severity"].astype(float).values

    # jitter in x to reduce overplotting
    rng = np.random.default_rng(42)
    xj = x + rng.normal(0.0, 0.06, size=len(x))

    fig, ax = plt.subplots(figsize=(10, 5))
    for r in regimes:
        m = df["primary_regime"] == r
        ax.scatter(
            xj[m.values],
            y[m.values],
            s=18,
            alpha=0.85,
            label=r,
        )

    ax.set_xticks(list(x_map.values()), list(x_map.keys()), rotation=0)
    ax.set_ylabel("maintenance_severity (objective / max objective per run)")
    ax.set_xlabel("primary_regime")
    ax.set_title("Regime map (colored by primary_regime; no run_id legend)")
    ax.grid(True, axis="y", alpha=0.25)
    ax.legend(title="primary_regime", fontsize=8, title_fontsize=9, loc="best", frameon=True)
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _resolve_anchor_features_parquet(results_root: Path, override_path: str | None) -> Path:
    if override_path:
        p = Path(override_path)
        if p.exists():
            return p
        raise FileNotFoundError(f"--anchor-features not found: {p}")

    # Preferred baseline run (as provided in the request)
    preferred = results_root / "campaigns" / "C3_o2_mid" / "run__acFree__o2lb-50__nh4cap-100" / "features.parquet"
    if preferred.exists():
        return preferred

    # Fallback: pick the first available features.parquet under campaigns/
    campaigns = results_root / "campaigns"
    if campaigns.exists():
        for p in campaigns.rglob("features.parquet"):
            return p

    raise FileNotFoundError(f"No features.parquet found under: {campaigns}")


def _resolve_calibrated_features_parquet(results_root: Path, override_path: str | None) -> Path:
    if override_path:
        p = Path(override_path)
        if p.exists():
            return p
        raise FileNotFoundError(f"--after-features not found: {p}")

    preferred = results_root / "campaigns" / "C6_atpm_calibrated" / "run__atpmCalib_linear" / "features.parquet"
    if preferred.exists():
        return preferred

    raise FileNotFoundError(
        "Calibrated run features.parquet not found. "
        "Expected: results/campaigns/C6_atpm_calibrated/run__atpmCalib_linear/features.parquet "
        "(under extracted results_root). You can also pass --after-features explicitly."
    )


def _plot_experiment_anchors_4panel(
    *,
    conditions_csv: Path,
    out_png: Path,
    out_csv: Path,
    figsize: tuple[int, int] = (10, 7),
) -> None:
    """
    Experimental anchors 4-panel (2x2) from conditions_experiment.csv.

    Outputs:
      - out_png: Fig02A_experiment_anchors_4panel.png
      - out_csv: Fig02A_experiment_anchors_4panel.csv (summarized points: mean/std/n)
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_png.parent.mkdir(parents=True, exist_ok=True)
    out_csv.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(conditions_csv)
    if "measured_OD" not in df.columns:
        raise ValueError("conditions_experiment.csv must include measured_OD")

    # Coerce numerics (robust to strings)
    for c in ["yeast_extract_gL", "pH0", "nh4cl_gL", "acetate_mM", "measured_OD"]:
        if c in df.columns:
            df[c] = pd.to_numeric(df[c], errors="coerce")

    if "set_name" in df.columns:
        df["set_name"] = df["set_name"].astype(str)

    def _subset(preferred_set: str, required: list[str]) -> pd.DataFrame:
        dd = df.copy()
        if "set_name" in dd.columns:
            cand = dd[dd["set_name"] == preferred_set].copy()
            if len(cand) > 0:
                dd = cand
        return dd.dropna(subset=required).copy()

    def _summarize(dd: pd.DataFrame, *, panel: str, x_col: str, group_col: str | None = None) -> pd.DataFrame:
        by = [x_col] + ([group_col] if group_col else [])
        agg = (
            dd.groupby(by, dropna=True)["measured_OD"]
            .agg(n="count", mean_OD="mean", std_OD="std")
            .reset_index()
            .sort_values(by=by)
        )
        agg.insert(0, "panel", panel)
        return agg

    yeast_dd = _subset("yeast_gradient", ["yeast_extract_gL", "measured_OD"])
    ph_dd = _subset("ph_toggle", ["pH0", "yeast_extract_gL", "measured_OD"])
    nh4_dd = _subset("nh4_gradient", ["nh4cl_gL", "measured_OD"])
    ac_dd = _subset("acetate_gradient", ["acetate_mM", "measured_OD"])

    rows: list[pd.DataFrame] = []
    rows.append(_summarize(yeast_dd, panel="Yeast gradient", x_col="yeast_extract_gL"))

    if len(ph_dd) > 0:
        ph_dd2 = ph_dd.copy()
        # Group yeast into absent/present bins; prefer exact 0.0 / 0.5 when available.
        ph_dd2["yeast_group_gL"] = np.where(ph_dd2["yeast_extract_gL"] >= 0.25, 0.5, 0.0)
        rows.append(_summarize(ph_dd2, panel="pH effect (yeast on/off)", x_col="pH0", group_col="yeast_group_gL"))

    rows.append(_summarize(nh4_dd, panel="NH4Cl gradient", x_col="nh4cl_gL"))
    rows.append(_summarize(ac_dd, panel="Acetate gradient", x_col="acetate_mM"))

    pd.concat(rows, ignore_index=True).to_csv(out_csv, index=False)

    fig, axes = plt.subplots(2, 2, figsize=figsize)
    fig.suptitle("Experimental anchors (OD600 at 32 h; n≈3)", fontsize=14)

    def _plot_line_marker(ax, dd: pd.DataFrame, *, x: str, title: str, xlabel: str, note: str | None = None) -> None:
        ax.set_title(title, fontsize=12)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("Measured OD600 at 32 h")
        ax.grid(True, alpha=0.25)
        if len(dd) == 0:
            ax.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax.transAxes)
            if note:
                ax.text(0.02, 0.02, note, ha="left", va="bottom", transform=ax.transAxes, fontsize=9)
            return
        g = dd.groupby(x, dropna=True)["measured_OD"].mean().reset_index().sort_values(x)
        ax.plot(g[x], g["measured_OD"], marker="o", lw=1.8)
        if note:
            ax.text(0.02, 0.02, note, ha="left", va="bottom", transform=ax.transAxes, fontsize=9)

    ax1, ax2, ax3, ax4 = axes[0, 0], axes[0, 1], axes[1, 0], axes[1, 1]

    _plot_line_marker(
        ax1,
        yeast_dd,
        x="yeast_extract_gL",
        title="Yeast gradient (set: yeast_gradient)",
        xlabel="Yeast extract (g/L)",
    )

    # pH effect with/without yeast
    ax2.set_title("pH effect with/without yeast (set: ph_toggle)", fontsize=12)
    ax2.set_xlabel("Initial pH (pH0)")
    ax2.set_ylabel("Measured OD600 at 32 h")
    ax2.grid(True, alpha=0.25)
    if len(ph_dd) == 0:
        ax2.text(0.5, 0.5, "No data", ha="center", va="center", transform=ax2.transAxes)
    else:
        ph_dd2 = ph_dd.copy()
        ph_dd2["yeast_group_gL"] = np.where(ph_dd2["yeast_extract_gL"] >= 0.25, 0.5, 0.0)
        for yv, label in [(0.5, "Yeast present(0.5 g/L)"), (0.0, "Yeast absent(0 g/L)")]:
            dd = ph_dd2[ph_dd2["yeast_group_gL"] == yv].copy()
            if len(dd) == 0:
                continue
            g = dd.groupby("pH0", dropna=True)["measured_OD"].mean().reset_index().sort_values("pH0")
            ax2.plot(g["pH0"], g["measured_OD"], marker="o", lw=1.8, label=label)
        ax2.legend(loc="best", fontsize=9, frameon=True)

    _plot_line_marker(
        ax3,
        nh4_dd,
        x="nh4cl_gL",
        title="NH4Cl gradient (set: nh4_gradient)",
        xlabel="NH4Cl (g/L)",
        note="YE fixed at 0.2 g/L",
    )
    _plot_line_marker(
        ax4,
        ac_dd,
        x="acetate_mM",
        title="Acetate gradient (set: acetate_gradient)",
        xlabel="Acetate (mM)",
    )

    fig.tight_layout(rect=[0, 0, 1, 0.95])
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)


def _plot_anchor_scatter(
    *,
    features_parquet: Path,
    out_png: Path,
    out_csv: Path | None,
    color_by_set: bool = False,
    figsize: tuple[int, int] = (6, 4),
) -> float:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from scipy.stats import spearmanr
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Missing dependency: scipy ({e}). Install: pip install scipy") from e

    out_png.parent.mkdir(parents=True, exist_ok=True)
    if out_csv is not None:
        out_csv.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_parquet(features_parquet)
    need = {"condition_id", "measured_OD", "set_name"}
    missing = sorted(need - set(df.columns))
    if missing:
        raise ValueError(f"features.parquet missing required columns: {missing}")

    d = df[list(need)].copy()

    # objective_value may not be present in features.parquet (wide feature matrix).
    # If missing, merge from sibling regime_fba.parquet in the same run folder.
    if "objective_value" in df.columns:
        d["objective_value"] = pd.to_numeric(df["objective_value"], errors="coerce")
    else:
        regime_fba = features_parquet.parent / "regime_fba.parquet"
        if not regime_fba.exists():
            raise ValueError(
                "features.parquet does not contain objective_value, and sibling regime_fba.parquet not found at: "
                f"{regime_fba}"
            )
        rf = pd.read_parquet(regime_fba, columns=["condition_id", "objective_value"])
        rf["objective_value"] = pd.to_numeric(rf["objective_value"], errors="coerce")
        d = d.merge(rf, on="condition_id", how="left")

    d["measured_OD"] = pd.to_numeric(d["measured_OD"], errors="coerce")
    d["set_name"] = d["set_name"].astype(str)
    d = d.dropna(subset=["objective_value", "measured_OD"]).copy()

    rho = float("nan")
    if len(d) >= 2:
        rho, _p = spearmanr(d["objective_value"].to_numpy(), d["measured_OD"].to_numpy())
        rho = float(rho)

    # Save data used for plotting
    if out_csv is not None:
        d[["condition_id", "objective_value", "measured_OD", "set_name"]].to_csv(out_csv, index=False)

    fig, ax = plt.subplots(figsize=figsize)
    if color_by_set:
        sets = sorted(d["set_name"].unique().tolist())
        for s in sets:
            dd = d[d["set_name"] == s]
            ax.scatter(dd["objective_value"], dd["measured_OD"], s=28, alpha=0.85, label=s)
        ax.legend(title="set_name", fontsize=8, title_fontsize=9, loc="best", frameon=True)
    else:
        ax.scatter(d["objective_value"], d["measured_OD"], s=28, alpha=0.85)

    ax.set_xlabel("Predicted growth (FBA objective)")
    ax.set_ylabel("Measured OD600 at 32 h")
    ax.grid(True, alpha=0.25)
    ax.text(0.98, 0.98, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, ha="right", va="top")
    fig.tight_layout()
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return rho


def _load_anchor_scatter_data(features_parquet: Path) -> pd.DataFrame:
    """
    Load the exact data used for the objective_value vs measured_OD scatter.
    Ensures objective_value exists by merging sibling regime_fba.parquet when needed.
    """
    df = pd.read_parquet(features_parquet)
    need = {"condition_id", "measured_OD", "set_name"}
    missing = sorted(need - set(df.columns))
    if missing:
        raise ValueError(f"features.parquet missing required columns: {missing}")

    d = df[list(need)].copy()
    if "objective_value" in df.columns:
        d["objective_value"] = pd.to_numeric(df["objective_value"], errors="coerce")
    else:
        regime_fba = features_parquet.parent / "regime_fba.parquet"
        if not regime_fba.exists():
            raise ValueError(
                "features.parquet does not contain objective_value, and sibling regime_fba.parquet not found at: "
                f"{regime_fba}"
            )
        rf = pd.read_parquet(regime_fba, columns=["condition_id", "objective_value"])
        rf["objective_value"] = pd.to_numeric(rf["objective_value"], errors="coerce")
        d = d.merge(rf, on="condition_id", how="left")

    d["measured_OD"] = pd.to_numeric(d["measured_OD"], errors="coerce")
    d["set_name"] = d["set_name"].astype(str)
    d["condition_id"] = d["condition_id"].astype(str)
    d = d.dropna(subset=["objective_value", "measured_OD"]).copy()
    return d[["condition_id", "objective_value", "measured_OD", "set_name"]]


def _plot_scatter_from_df(
    df: pd.DataFrame,
    *,
    out_png: Path,
    title: str | None = None,
    figsize: tuple[int, int] = (6, 4),
) -> float:
    """
    Plot objective_value vs measured_OD from a prepared dataframe.
    Required columns: condition_id, objective_value, measured_OD, set_name.
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from scipy.stats import spearmanr
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Missing dependency: scipy ({e}). Install: pip install scipy") from e

    if df.empty:
        raise ValueError("No rows to plot (df is empty).")

    d = df.copy()
    d["objective_value"] = pd.to_numeric(d["objective_value"], errors="coerce")
    d["measured_OD"] = pd.to_numeric(d["measured_OD"], errors="coerce")
    d = d.dropna(subset=["objective_value", "measured_OD"]).copy()
    if len(d) < 2:
        raise ValueError("Need at least 2 valid rows to compute Spearman rho.")

    rho, _p = spearmanr(d["objective_value"].to_numpy(), d["measured_OD"].to_numpy())
    rho = float(rho)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(d["objective_value"], d["measured_OD"], s=28, alpha=0.85)
    ax.set_xlabel("Predicted growth (FBA objective)")
    ax.set_ylabel("Measured OD600 at 32 h")
    if title:
        ax.set_title(title)
    ax.grid(True, alpha=0.25)
    ax.text(0.98, 0.98, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, ha="right", va="top")
    fig.tight_layout()
    out_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close(fig)
    return rho


def _load_feature_matrix_regime_dataset(regime_dataset_parquet: Path, *, seed: int = 42, max_rows: int = 2000) -> pd.DataFrame:
    """
    Load X feature matrix from regime_dataset.parquet:
      - select width__/mid__/signchange__ columns
      - coerce to numeric
      - sample up to max_rows for faster SHAP plotting
    """
    df = pd.read_parquet(regime_dataset_parquet)
    feat_cols = [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]
    if not feat_cols:
        raise ValueError("No feature columns found in regime_dataset.parquet (expected width__/mid__/signchange__).")
    X = df[feat_cols].copy()
    for c in X.columns:
        if X[c].dtype == bool:
            X[c] = X[c].astype(int)
    X = X.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    if len(X) > max_rows:
        X = X.sample(n=int(max_rows), random_state=int(seed))
    return X


def _plot_regime_vs_od_and_severity(
    *,
    regime_dataset_parquet: Path,
    out_od_png: Path,
    out_od_csv: Path,
    out_severity_png: Path | None = None,
) -> None:
    """
    Fig02B/Fig02C:
      - Fig02B: boxplot + jitter for measured_OD by primary_regime
      - Fig02C (optional): boxplot + jitter for maintenance_severity by primary_regime

    Input: regime_dataset.parquet (preferred)
      expected cols: primary_regime, measured_OD, maintenance_severity (optional), run_id (optional), condition_id (optional)
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rng = np.random.default_rng(42)

    df = pd.read_parquet(regime_dataset_parquet)

    # Column normalization: the dataset may have measured_OD_x / measured_OD_y after merges.
    def _pick_col(candidates: list[str]) -> str | None:
        for c in candidates:
            if c in df.columns:
                return c
        return None

    regime_col = _pick_col(["primary_regime", "label"])
    od_col = _pick_col(["measured_OD", "measured_OD_y", "measured_OD_x"])
    sev_col = _pick_col(["maintenance_severity"])
    run_col = _pick_col(["run_id"])
    cond_col = _pick_col(["condition_id"])

    if regime_col is None:
        raise ValueError("regime_dataset.parquet must include one of: primary_regime, label")
    if od_col is None:
        raise ValueError("regime_dataset.parquet must include one of: measured_OD, measured_OD_y, measured_OD_x")

    d = df.copy()
    d["primary_regime"] = d[regime_col].astype(str).str.strip()
    d["measured_OD"] = pd.to_numeric(d[od_col], errors="coerce")

    # Compute maintenance_severity if missing but objective_value + run_id exist.
    if sev_col is not None:
        d["maintenance_severity"] = pd.to_numeric(d[sev_col], errors="coerce")
    else:
        if "objective_value" in d.columns and run_col is not None:
            obj = pd.to_numeric(d["objective_value"], errors="coerce")
            d["_obj"] = obj
            max_by_run = d.groupby(run_col, dropna=False)["_obj"].transform("max")
            d["maintenance_severity"] = d["_obj"] / max_by_run
            d.loc[~np.isfinite(d["maintenance_severity"]), "maintenance_severity"] = np.nan
            d.drop(columns=["_obj"], inplace=True)
        else:
            d["maintenance_severity"] = np.nan

    # Keep only rows with OD for Fig02B
    d_od = d.dropna(subset=["primary_regime", "measured_OD"]).copy()

    # Category order (keep empty categories)
    regimes_order = ["Ac_limited", "N_limited", "Pi_limited", "O2_limited", "Unconstrained"]
    # Keep additional regimes if present
    extra = [r for r in sorted(d_od["primary_regime"].unique().tolist()) if r not in regimes_order]
    regimes = regimes_order + extra

    # Save data used for plotting
    cols = ["primary_regime", "measured_OD", "maintenance_severity"]
    if run_col is not None:
        d_od["run_id"] = d_od[run_col].astype(str)
        cols.append("run_id")
    if cond_col is not None:
        d_od["condition_id"] = d_od[cond_col].astype(str)
        cols.append("condition_id")
    out_od_png.parent.mkdir(parents=True, exist_ok=True)
    out_od_csv.parent.mkdir(parents=True, exist_ok=True)
    d_od[cols].to_csv(out_od_csv, index=False)

    # Simple color map (no legend spam)
    palette = {
        "Ac_limited": "#d62728",
        "O2_limited": "#1f77b4",
        "N_limited": "#2ca02c",
        "Pi_limited": "#9467bd",
        "Unconstrained": "#7f7f7f",
    }

    def _draw_box_jitter(ax, y_col: str, title: str, ylabel: str) -> None:
        ax.set_title(title, fontsize=12)
        ax.set_ylabel(ylabel)
        ax.grid(True, axis="y", alpha=0.25)
        positions = np.arange(1, len(regimes) + 1)

        data = []
        ns = []
        for r in regimes:
            vals = d_od.loc[d_od["primary_regime"] == r, y_col].dropna().to_numpy()
            data.append(vals)
            ns.append(int(len(vals)))

        # Matplotlib boxplot tolerates empty arrays; showfliers=False to reduce clutter
        bp = ax.boxplot(
            data,
            positions=positions,
            widths=0.6,
            patch_artist=True,
            showfliers=False,
            medianprops=dict(color="black", linewidth=1.2),
        )
        for patch, r in zip(bp["boxes"], regimes):
            patch.set_facecolor(palette.get(r, "#bbbbbb"))
            patch.set_alpha(0.35)
            patch.set_edgecolor("#444444")

        # Jittered points
        for i, (r, vals) in enumerate(zip(regimes, data), start=1):
            if len(vals) == 0:
                continue
            xj = i + rng.normal(0, 0.06, size=len(vals))
            ax.scatter(
                xj,
                vals,
                s=18,
                alpha=0.7,
                color=palette.get(r, "#555555"),
                edgecolors="none",
            )

        ax.set_xticks(positions)
        ax.set_xticklabels(regimes, rotation=0)

        # n annotations
        finite = d_od[y_col].dropna()
        if len(finite) > 0:
            y_max = float(finite.max())
            y_min = float(finite.min())
            y_span = max(1e-9, y_max - y_min)
            y_text = y_max + 0.06 * y_span
            for x, n in zip(positions, ns):
                ax.text(x, y_text, f"n={n}", ha="center", va="bottom", fontsize=9, color="#333333")
            ax.set_ylim(y_min - 0.05 * y_span, y_max + 0.14 * y_span)

    # Fig02B: OD vs regime
    fig, ax = plt.subplots(figsize=(10, 4.8))
    _draw_box_jitter(
        ax,
        y_col="measured_OD",
        title="Experimental OD vs predicted limiting regime",
        ylabel="Measured OD600 at 32 h",
    )
    fig.tight_layout()
    fig.savefig(out_od_png, dpi=220, bbox_inches="tight")
    plt.close(fig)

    # Fig02C: severity vs regime (optional)
    if out_severity_png is not None and "maintenance_severity" in d.columns:
        d_sev = d.dropna(subset=["primary_regime", "maintenance_severity"]).copy()
        if len(d_sev) > 0:
            # Use same regime list but severity data subset
            d_od_backup = d_od
            d_od = d_sev  # reuse helper
            fig2, ax2 = plt.subplots(figsize=(10, 4.2))
            _draw_box_jitter(
                ax2,
                y_col="maintenance_severity",
                title="Maintenance severity vs predicted limiting regime",
                ylabel="Maintenance severity (objective / run max)",
            )
            fig2.tight_layout()
            out_severity_png.parent.mkdir(parents=True, exist_ok=True)
            fig2.savefig(out_severity_png, dpi=220, bbox_inches="tight")
            plt.close(fig2)
            d_od = d_od_backup


def _save_group_summary(
    df: pd.DataFrame,
    *,
    group_col: str,
    y_col: str,
    out_csv: Path,
    regimes: list[str],
) -> None:
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for r in regimes:
        vals = pd.to_numeric(df.loc[df[group_col] == r, y_col], errors="coerce").dropna()
        if len(vals) == 0:
            continue
        q25 = float(vals.quantile(0.25))
        q75 = float(vals.quantile(0.75))
        rows.append(
            {
                "primary_regime": r,
                "n_total": int(len(vals)),
                "median": float(vals.median()),
                "q25": q25,
                "q75": q75,
                "iqr": float(q75 - q25),
                "mean": float(vals.mean()),
                "std": float(vals.std(ddof=1)) if len(vals) >= 2 else float("nan"),
            }
        )
    pd.DataFrame(rows).to_csv(out_csv, index=False)


def _plot_regime_box_jitter_v2(
    df: pd.DataFrame,
    *,
    y_col: str,
    out_png: Path,
    out_summary_csv: Path,
    title: str,
    ylabel: str,
    xlabel: str,
    regimes: list[str],
    caption_not_observed: str | None,
    o2_subsample_n: int = 50,
    seed: int = 42,
    figsize: tuple[int, int] = (10, 5),
) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_png.parent.mkdir(parents=True, exist_ok=True)

    d = df.copy()
    d["primary_regime"] = d["primary_regime"].astype(str).str.strip()
    d[y_col] = pd.to_numeric(d[y_col], errors="coerce")
    d = d.dropna(subset=["primary_regime", y_col]).copy()

    # restrict to requested regimes, and drop empty regimes from plot
    d = d[d["primary_regime"].isin(regimes)].copy()
    present = [r for r in regimes if (d["primary_regime"] == r).any()]

    # Save summary for present regimes only
    _save_group_summary(d, group_col="primary_regime", y_col=y_col, out_csv=out_summary_csv, regimes=present)

    # Colors (no legend)
    palette = {
        "Ac_limited": "#d62728",
        "O2_limited": "#1f77b4",
        "N_limited": "#2ca02c",
    }

    rng = np.random.default_rng(int(seed))
    positions = np.arange(1, len(present) + 1)

    data = []
    ns = []
    for r in present:
        vals = d.loc[d["primary_regime"] == r, y_col].dropna().to_numpy()
        data.append(vals)
        ns.append(int(len(vals)))

    fig, ax = plt.subplots(figsize=figsize)
    ax.set_title(title, fontsize=13)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y", alpha=0.22)

    bp = ax.boxplot(
        data,
        positions=positions,
        widths=0.65,
        patch_artist=True,
        showfliers=False,
        medianprops=dict(color="black", linewidth=1.3),
    )
    for patch, r in zip(bp["boxes"], present):
        patch.set_facecolor(palette.get(r, "#bbbbbb"))
        patch.set_alpha(0.28)
        patch.set_edgecolor("#444444")

    # scatter overlay (subsample only for O2_limited)
    for i, (r, vals) in enumerate(zip(present, data), start=1):
        if len(vals) == 0:
            continue
        vals_plot = vals
        if r == "O2_limited" and len(vals) > int(o2_subsample_n):
            idx = rng.choice(len(vals), size=int(o2_subsample_n), replace=False)
            vals_plot = vals[idx]
        xj = i + rng.normal(0, 0.07, size=len(vals_plot))
        ax.scatter(
            xj,
            vals_plot,
            s=20,
            alpha=0.35,
            color=palette.get(r, "#555555"),
            edgecolors="none",
        )

    ax.set_xticks(positions)
    ax.set_xticklabels(present)

    # n labels: show total n (not subsample)
    finite = pd.to_numeric(d[y_col], errors="coerce").dropna()
    if len(finite) > 0:
        y_max = float(finite.max())
        y_min = float(finite.min())
        y_span = max(1e-9, y_max - y_min)
        y_text = y_max + 0.06 * y_span
        for x, n in zip(positions, ns):
            ax.text(x, y_text, f"n={n}", ha="center", va="bottom", fontsize=10, color="#333333")
        ax.set_ylim(y_min - 0.05 * y_span, y_max + 0.16 * y_span)

    if caption_not_observed:
        fig.text(0.01, 0.01, caption_not_observed, ha="left", va="bottom", fontsize=10, color="#333333")

    fig.tight_layout(rect=[0, 0.03, 1, 1])
    fig.savefig(out_png, dpi=260, bbox_inches="tight")
    plt.close(fig)


def _load_regime_dataset_normalized(regime_dataset_parquet: Path) -> pd.DataFrame:
    """
    Normalize regime_dataset.parquet into a consistent schema:
      - primary_regime (from primary_regime or label)
      - measured_OD (from measured_OD or measured_OD_y/x)
      - maintenance_severity (use existing if present; else compute from objective_value/run_id when possible)
    Keeps feature columns (width__/mid__/signchange__) if present.
    """
    df = pd.read_parquet(regime_dataset_parquet)

    def _pick_col(_df: pd.DataFrame, candidates: list[str]) -> str | None:
        for c in candidates:
            if c in _df.columns:
                return c
        return None

    regime_col = _pick_col(df, ["primary_regime", "label"])
    od_col = _pick_col(df, ["measured_OD", "measured_OD_y", "measured_OD_x"])
    run_col = _pick_col(df, ["run_id"])
    cond_col = _pick_col(df, ["condition_id"])

    if regime_col is None:
        raise ValueError("regime_dataset.parquet missing regime column (expected primary_regime or label)")
    if od_col is None:
        raise ValueError("regime_dataset.parquet missing OD column (expected measured_OD or measured_OD_y/x)")

    out = df.copy()
    out["primary_regime"] = out[regime_col].astype(str).str.strip()
    out["measured_OD"] = pd.to_numeric(out[od_col], errors="coerce")

    if "maintenance_severity" in out.columns:
        out["maintenance_severity"] = pd.to_numeric(out["maintenance_severity"], errors="coerce")
    else:
        if "objective_value" in out.columns and run_col is not None:
            obj = pd.to_numeric(out["objective_value"], errors="coerce")
            out["_obj"] = obj
            max_by_run = out.groupby(run_col, dropna=False)["_obj"].transform("max")
            out["maintenance_severity"] = out["_obj"] / max_by_run
            out.loc[~np.isfinite(out["maintenance_severity"]), "maintenance_severity"] = np.nan
            out.drop(columns=["_obj"], inplace=True)
        else:
            out["maintenance_severity"] = np.nan

    if run_col is not None and run_col != "run_id":
        out["run_id"] = out[run_col].astype(str)
    if cond_col is not None and cond_col != "condition_id":
        out["condition_id"] = out[cond_col].astype(str)

    return out


def _plot_od_vs_severity_scatter(
    df: pd.DataFrame,
    *,
    out_png: Path,
    figsize: tuple[int, int] = (6, 4),
) -> float:
    """
    Fig02D: measured_OD vs maintenance_severity scatter + Spearman rho.
    Returns rho (float).
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    try:
        from scipy.stats import spearmanr
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Missing dependency: scipy ({e}). Install: pip install scipy") from e

    out_png.parent.mkdir(parents=True, exist_ok=True)

    d = df.copy()
    d["measured_OD"] = pd.to_numeric(d["measured_OD"], errors="coerce")
    d["maintenance_severity"] = pd.to_numeric(d["maintenance_severity"], errors="coerce")
    d = d.dropna(subset=["measured_OD", "maintenance_severity"]).copy()
    if len(d) < 2:
        raise ValueError("Not enough rows with measured_OD and maintenance_severity to plot Fig02D.")

    rho, _p = spearmanr(d["maintenance_severity"].to_numpy(), d["measured_OD"].to_numpy())
    rho = float(rho)

    fig, ax = plt.subplots(figsize=figsize)
    ax.scatter(d["maintenance_severity"], d["measured_OD"], s=26, alpha=0.55, edgecolors="none")
    ax.set_xlabel("Maintenance severity (objective / run max)")
    ax.set_ylabel("Measured OD600 at 32 h")
    ax.set_title("Measured OD vs maintenance severity", fontsize=12)
    ax.grid(True, alpha=0.25)
    ax.text(0.98, 0.98, f"Spearman ρ = {rho:.2f}", transform=ax.transAxes, ha="right", va="top")
    fig.tight_layout()
    fig.savefig(out_png, dpi=260, bbox_inches="tight")
    plt.close(fig)
    return rho


def _train_od_low_classifier_and_report(
    df: pd.DataFrame,
    *,
    out_shap_png: Path,
    out_metrics_csv: Path,
    out_report_csv: Path,
    od_threshold: float = 0.6,
    seed: int = 42,
) -> None:
    """
    OD_low classifier (OD < threshold) with a SHAP-like explanation plot.

    Practical note: We compute "SHAP-like" per-feature contributions for a linear model:
      contrib = (x - mean_train) * coef
    and plot mean(|contrib|) across test samples (bar plot).
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import accuracy_score, roc_auc_score
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    out_shap_png.parent.mkdir(parents=True, exist_ok=True)
    out_metrics_csv.parent.mkdir(parents=True, exist_ok=True)
    out_report_csv.parent.mkdir(parents=True, exist_ok=True)

    # Aggregate to condition-level (22 conditions): stabilize across multiple runs
    feat_cols = [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]
    base_cols = ["condition_id", "measured_OD", "primary_regime", "maintenance_severity"]
    have = [c for c in base_cols if c in df.columns]

    d0 = df[have + feat_cols].copy()
    d0["measured_OD"] = pd.to_numeric(d0["measured_OD"], errors="coerce")
    d0["maintenance_severity"] = pd.to_numeric(d0["maintenance_severity"], errors="coerce")

    # groupby condition_id: OD is constant; severity/feats use median/mean
    agg = {"measured_OD": "mean", "maintenance_severity": "median", "primary_regime": lambda s: s.dropna().astype(str).mode().iloc[0] if len(s.dropna()) else ""}
    for c in feat_cols:
        agg[c] = "mean"
    d = d0.groupby("condition_id", dropna=False).agg(agg).reset_index()

    # Label
    d["OD_low"] = (d["measured_OD"] < float(od_threshold)).astype(int)

    # Pick a small feature set: severity + top correlated features with OD_low
    cand_cols = ["maintenance_severity"] + feat_cols
    X_all = d[cand_cols].copy()
    for c in X_all.columns:
        X_all[c] = pd.to_numeric(X_all[c], errors="coerce")
    X_all = X_all.fillna(0.0)
    y_all = d["OD_low"].astype(int).to_numpy()

    # Correlation-based selection (simple + stable)
    corrs = []
    y_center = y_all - float(np.mean(y_all))
    for c in X_all.columns:
        x = X_all[c].to_numpy()
        x_center = x - float(np.mean(x))
        denom = float(np.sqrt(np.sum(x_center**2) * np.sum(y_center**2)))
        corr = float(np.sum(x_center * y_center) / denom) if denom > 0 else 0.0
        corrs.append((abs(corr), c))
    corrs.sort(reverse=True)
    topk = 8  # keep small for interpretability
    selected = []
    for _a, c in corrs:
        if c == "maintenance_severity":
            continue
        selected.append(c)
        if len(selected) >= topk:
            break
    feature_cols = ["maintenance_severity"] + selected

    X = X_all[feature_cols].copy()

    # Split by conditions (22 rows): stratified if possible
    test_size = 0.3 if len(d) >= 10 else 0.2
    X_train, X_test, y_train, y_test, idx_train, idx_test = train_test_split(
        X,
        y_all,
        np.arange(len(d)),
        test_size=float(test_size),
        random_state=int(seed),
        stratify=y_all if len(np.unique(y_all)) > 1 else None,
    )

    pipe = Pipeline(
        steps=[
            ("scaler", StandardScaler()),
            ("clf", LogisticRegression(max_iter=5000, random_state=int(seed))),
        ]
    )
    pipe.fit(X_train, y_train)
    proba = pipe.predict_proba(X_test)[:, 1]
    pred = (proba >= 0.5).astype(int)

    auc = float("nan")
    if len(np.unique(y_test)) > 1:
        auc = float(roc_auc_score(y_test, proba))
    acc = float(accuracy_score(y_test, pred))

    pd.DataFrame(
        [
            {
                "auc": auc,
                "acc": acc,
                "n_conditions": int(len(d)),
                "n_train": int(len(y_train)),
                "n_test": int(len(y_test)),
                "od_threshold": float(od_threshold),
                "features_used": ";".join(feature_cols),
            }
        ]
    ).to_csv(out_metrics_csv, index=False)

    # SHAP-like contributions on test set: (x - mean_train) * coef (in standardized space -> use pipeline pieces)
    scaler: StandardScaler = pipe.named_steps["scaler"]
    clf: LogisticRegression = pipe.named_steps["clf"]
    X_train_z = scaler.transform(X_train)
    X_test_z = scaler.transform(X_test)

    coef = clf.coef_.reshape(-1)  # (F,)
    mu_train = np.mean(X_train_z, axis=0)
    contrib = (X_test_z - mu_train) * coef  # (N,F)
    mean_abs = np.mean(np.abs(contrib), axis=0)

    order = np.argsort(mean_abs)[::-1]
    topn = min(10, len(feature_cols))
    top_idx = order[:topn]
    top_feats = [feature_cols[i] for i in top_idx]
    top_vals = mean_abs[top_idx]

    # Plot bar (paper-friendly, no legend explosion)
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(9, 4.8))
    y_pos = np.arange(len(top_feats))[::-1]
    ax.barh(y_pos, top_vals[::-1], color="#4c72b0", alpha=0.85)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(top_feats[::-1], fontsize=9)
    ax.set_xlabel("Mean |contribution| (SHAP-like, logistic)")
    ax.set_title("OD_low classifier: top feature contributions", fontsize=12)
    ax.text(0.98, 0.02, f"AUC={auc:.2f}, Acc={acc:.2f}", transform=ax.transAxes, ha="right", va="bottom", fontsize=10)
    fig.tight_layout()
    fig.savefig(out_shap_png, dpi=260, bbox_inches="tight")
    plt.close(fig)

    # Build 22-condition diagnosis report (fit final model on ALL conditions for stable per-condition contributions)
    pipe.fit(X, y_all)
    scaler2: StandardScaler = pipe.named_steps["scaler"]
    clf2: LogisticRegression = pipe.named_steps["clf"]
    # Keep feature names (avoid sklearn warning)
    X_z = scaler2.transform(X)
    coef2 = clf2.coef_.reshape(-1)
    mu2 = np.mean(X_z, axis=0)
    contrib_all = (X_z - mu2) * coef2  # (22,F)

    sev = pd.to_numeric(d["maintenance_severity"], errors="coerce")
    sev_pct = sev.rank(pct=True) * 100.0

    top5_list = []
    for i in range(len(d)):
        row = contrib_all[i, :]
        idx = np.argsort(np.abs(row))[::-1][:5]
        items = [f"{feature_cols[j]}={row[j]:.3g}" for j in idx]
        top5_list.append(";".join(items))

    report = pd.DataFrame(
        {
            "condition_id": d["condition_id"].astype(str),
            "measured_OD": d["measured_OD"].astype(float),
            "primary_regime": d["primary_regime"].astype(str),
            "maintenance_severity": sev.astype(float),
            "severity_percentile": sev_pct.astype(float),
            "top5_shap_features": top5_list,
        }
    )
    report.to_csv(out_report_csv, index=False)


def _plot_shap_beeswarm_from_saved_model(
    *,
    model_json: Path,
    X: pd.DataFrame,
    out_png: Path,
    is_classifier: bool,
    max_display: int = 15,
    figsize: tuple[int, int] = (10, 5),
) -> None:
    """
    Recompute SHAP beeswarm from saved xgboost model.json and feature matrix X.
    """
    try:
        import shap
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Missing dependency: shap ({e}). Install: pip install shap") from e
    try:
        from xgboost import XGBClassifier, XGBRegressor
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Missing dependency: xgboost ({e}). Install: pip install xgboost") from e

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_png.parent.mkdir(parents=True, exist_ok=True)

    if is_classifier:
        model = XGBClassifier()
    else:
        model = XGBRegressor()
    model.load_model(str(model_json))

    # Prefer SHAP TreeExplainer, but fall back to xgboost pred_contribs when TreeExplainer fails
    # (e.g., some multiclass JSON models store base_score as a list string).
    try:
        explainer = shap.TreeExplainer(model)
        shap_values = explainer.shap_values(X)
    except Exception as e:  # noqa: BLE001
        # xgboost-native SHAP values (no TreeExplainer parsing)
        import xgboost as xgb

        print(f"[WARN] shap.TreeExplainer failed ({e}). Falling back to xgboost pred_contribs.")
        booster = model.get_booster()
        dmat = xgb.DMatrix(X, feature_names=list(X.columns))
        try:
            contrib = booster.predict(dmat, pred_contribs=True, strict_shape=True)
        except TypeError:
            contrib = booster.predict(dmat, pred_contribs=True)

        arr = np.asarray(contrib)
        # Regression: (N, F+1) -> drop bias
        if arr.ndim == 2:
            shap_values = arr[:, :-1]
        # Multiclass strict_shape: (N, C, F+1) -> list[C] of (N,F)
        elif arr.ndim == 3:
            shap_values = [arr[:, c, :-1] for c in range(arr.shape[1])]
        else:
            raise RuntimeError(f"Unexpected pred_contribs shape: {arr.shape}") from e

    plt.figure(figsize=figsize)
    shap.summary_plot(shap_values, X, show=False, max_display=int(max_display))
    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close()


def _rewrite_png_with_figsize(src_png: Path, out_png: Path, figsize: tuple[int, int] = (10, 5)) -> None:
    """
    Fallback when shap/xgboost are not available locally:
    re-save an existing PNG into a controlled canvas size.
    (Note: does not change max_display; only changes the figure size.)
    """
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.image as mpimg

    out_png.parent.mkdir(parents=True, exist_ok=True)
    img = mpimg.imread(src_png)
    plt.figure(figsize=figsize)
    plt.imshow(img)
    plt.axis("off")
    plt.tight_layout()
    plt.savefig(out_png, dpi=220, bbox_inches="tight")
    plt.close()


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Generate draft paper figures from extracted results.")
    parser.add_argument(
        "--conditions-csv",
        default=str(Path("acetate_xai") / "data" / "conditions_experiment.csv"),
        help="Experimental conditions CSV (must include measured_OD). Used for Fig02A experimental anchors panel.",
    )
    parser.add_argument(
        "--anchor-features",
        default=None,
        help="Path to baseline run features.parquet for Fig02E anchor scatter. "
        "If omitted, tries the default C3_o2_mid baseline run, then falls back to the first campaigns/**/features.parquet.",
    )
    parser.add_argument(
        "--before-features",
        default=None,
        help="Optional override for Fig02E BEFORE (objective vs OD). Defaults to the baseline C3_o2_mid run if available.",
    )
    parser.add_argument(
        "--after-features",
        default=None,
        help="Optional override for Fig02E AFTER (calibrated objective vs OD). Defaults to C6_atpm_calibrated run if available.",
    )
    parser.add_argument(
        "--anchor-color-by-set",
        action="store_true",
        help="Color Fig02E scatter points by set_name (adds legend).",
    )
    args = parser.parse_args(argv)

    results_root = _find_results_root()

    # New outputs only
    figures_out = Path("results") / "figures_draft"
    doc_out = Path("docs") / "PROJECT_STATUS_SUMMARY.md"

    # Source files inside results_root
    src_regime_map_png = results_root / "platform_summary" / "regime_map.png"
    src_regime_map_csv = results_root / "platform_summary" / "regime_map.csv"
    src_platform_csv = results_root / "platform_summary" / "platform_summary.csv"
    src_robustness_txt = results_root / "platform_summary" / "robustness_summary.txt"
    src_xgb_regime_model = results_root / "xai_xgb" / "model.json"
    src_xgb_sev_model = results_root / "xai_xgb_maintenance" / "model.json"
    src_shap_regime_png = results_root / "xai_xgb" / "shap_beeswarm.png"
    src_shap_sev_png = results_root / "xai_xgb_maintenance" / "shap_beeswarm.png"
    src_causal_png = results_root / "causal_dag" / "causal_dag.png"
    src_edge_stab = results_root / "causal_dag" / "edge_stability.csv"
    src_regime_dataset = results_root / "regime_dataset.parquet"

    required = [
        src_regime_map_png,
        src_regime_map_csv,
        src_platform_csv,
        src_robustness_txt,
        src_xgb_regime_model,
        src_xgb_sev_model,
        src_shap_regime_png,
        src_shap_sev_png,
        src_causal_png,
        src_edge_stab,
        src_regime_dataset,
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

    # Fig02A: experimental anchors (2x2 panels) from conditions_experiment.csv
    _plot_experiment_anchors_4panel(
        conditions_csv=Path(args.conditions_csv),
        out_png=figures_out / "Fig02A_experiment_anchors_4panel.png",
        out_csv=figures_out / "Fig02A_experiment_anchors_4panel.csv",
        figsize=(10, 7),
    )

    # Fig02B (and optional Fig02C): regime vs OD (and severity) from regime_dataset.parquet
    _plot_regime_vs_od_and_severity(
        regime_dataset_parquet=src_regime_dataset,
        out_od_png=figures_out / "Fig02B_regime_vs_OD.png",
        out_od_csv=figures_out / "Fig02B_regime_vs_OD.csv",
        out_severity_png=figures_out / "Fig02C_regime_vs_severity.png",
    )

    # Fig02B v2 / Fig02C v2 (paper-ready): 3 regimes only, subsample O2_limited scatter, add caption for non-observed.
    df_reg = pd.read_parquet(src_regime_dataset)
    # Normalize columns similar to _plot_regime_vs_od_and_severity
    def _pick_col(_df: pd.DataFrame, candidates: list[str]) -> str | None:
        for c in candidates:
            if c in _df.columns:
                return c
        return None

    _regime_col = _pick_col(df_reg, ["primary_regime", "label"])
    _od_col = _pick_col(df_reg, ["measured_OD", "measured_OD_y", "measured_OD_x"])
    _run_col = _pick_col(df_reg, ["run_id"])
    if _regime_col is None or _od_col is None:
        raise ValueError("regime_dataset.parquet missing regime/OD columns needed for Fig02B v2")

    dv2 = pd.DataFrame(
        {
            "primary_regime": df_reg[_regime_col].astype(str).str.strip(),
            "measured_OD": pd.to_numeric(df_reg[_od_col], errors="coerce"),
        }
    )
    if _run_col is not None and "objective_value" in df_reg.columns:
        obj = pd.to_numeric(df_reg["objective_value"], errors="coerce")
        tmp = pd.DataFrame({"run_id": df_reg[_run_col].astype(str), "objective_value": obj})
        max_by_run = tmp.groupby("run_id", dropna=False)["objective_value"].transform("max")
        dv2["maintenance_severity"] = obj / max_by_run
        dv2.loc[~np.isfinite(dv2["maintenance_severity"]), "maintenance_severity"] = np.nan
    else:
        dv2["maintenance_severity"] = pd.to_numeric(df_reg.get("maintenance_severity", np.nan), errors="coerce")

    regimes3 = ["Ac_limited", "N_limited", "O2_limited"]
    caption = "Pi_limited/Unconstrained not observed in current campaigns"

    _plot_regime_box_jitter_v2(
        dv2,
        y_col="measured_OD",
        out_png=figures_out / "Fig02B_regime_vs_OD_v2.png",
        out_summary_csv=figures_out / "Fig02B_regime_vs_OD_summary.csv",
        title="Experimental OD vs predicted limiting regime",
        ylabel="Measured OD600 at 32 h",
        xlabel="Primary limiting regime (FBA saturation label)",
        regimes=regimes3,
        caption_not_observed=caption,
        o2_subsample_n=50,
        seed=42,
        figsize=(10, 5),
    )

    _plot_regime_box_jitter_v2(
        dv2,
        y_col="maintenance_severity",
        out_png=figures_out / "Fig02C_regime_vs_severity.png",
        out_summary_csv=figures_out / "Fig02C_regime_vs_severity_summary.csv",
        title="Maintenance severity vs predicted limiting regime",
        ylabel="Maintenance severity (objective / run max)",
        xlabel="Primary limiting regime (FBA saturation label)",
        regimes=regimes3,
        caption_not_observed=caption,
        o2_subsample_n=50,
        seed=42,
        figsize=(10, 5),
    )

    # Fig02D + Fig02E + diagnosis report (diagnostic, not OD prediction)
    reg_norm = _load_regime_dataset_normalized(src_regime_dataset)
    _plot_od_vs_severity_scatter(
        reg_norm,
        out_png=figures_out / "Fig02D_OD_vs_severity.png",
        figsize=(6, 4),
    )
    _train_od_low_classifier_and_report(
        reg_norm,
        out_shap_png=figures_out / "Fig02E_low_growth_classifier_shap.png",
        out_metrics_csv=figures_out / "OD_low_metrics.csv",
        out_report_csv=figures_out / "diagnosis_report_22conds.csv",
        od_threshold=0.6,
        seed=42,
    )

    # Fig02: re-draw without run_id legend (primary_regime colors only)
    _plot_regime_map_primary_only(src_regime_map_csv, figures_out / "Fig02_regime_map.png")

    # Fig02E: anchor scatter (legacy single panel output)
    anchor_features = _resolve_anchor_features_parquet(results_root, args.anchor_features)
    _plot_anchor_scatter(
        features_parquet=anchor_features,
        out_png=figures_out / "Fig02E_anchor_scatter.png",
        out_csv=figures_out / "Fig02E_anchor_scatter_data.csv",
        color_by_set=bool(args.anchor_color_by_set),
        figsize=(6, 4),
    )

    # Fig02E before/after (objective vs OD; Spearman rho) for ATPM calibration evaluation
    before_features = _resolve_anchor_features_parquet(results_root, args.before_features)
    after_features = _resolve_calibrated_features_parquet(results_root, args.after_features)
    before_df = _load_anchor_scatter_data(before_features)
    after_df = _load_anchor_scatter_data(after_features)
    before_df.to_csv(figures_out / "Fig02E_before_objective_vs_OD_data.csv", index=False)
    after_df.to_csv(figures_out / "Fig02E_after_calibObjective_vs_OD_data.csv", index=False)

    # Also save a joined table (condition_id-aligned) for easy downstream analysis.
    joined = before_df.merge(
        after_df[["condition_id", "objective_value"]].rename(columns={"objective_value": "objective_value_after"}),
        on="condition_id",
        how="inner",
    ).rename(columns={"objective_value": "objective_value_before"})
    joined.to_csv(figures_out / "Fig02E_before_after_join.csv", index=False)

    rho_before = _plot_anchor_scatter(
        features_parquet=before_features,
        out_png=figures_out / "Fig02E_before_objective_vs_OD.png",
        out_csv=None,
        color_by_set=False,
        figsize=(6, 4),
    )
    rho_after = _plot_anchor_scatter(
        features_parquet=after_features,
        out_png=figures_out / "Fig02E_after_calibObjective_vs_OD.png",
        out_csv=None,
        color_by_set=False,
        figsize=(6, 4),
    )
    (figures_out / "Fig02E_rho_compare.txt").write_text(
        f"rho_before={rho_before:.6g}\n" f"rho_after={rho_after:.6g}\n",
        encoding="utf-8",
    )

    # Fig02E: ATPM calibration evaluation on acetate-gradient subset only (train/test split)
    # - subset: set_name == acetate_gradient
    # - train anchors: AC_25, AC_150
    # - test: AC_50, AC_100
    try:
        subset = after_df[after_df["set_name"] == "acetate_gradient"].copy()
        train_ids = {"AC_25", "AC_150"}
        test_ids = {"AC_50", "AC_100"}

        df_train = subset[subset["condition_id"].isin(train_ids)].copy()
        df_test = subset[subset["condition_id"].isin(test_ids)].copy()

        rho_train = _plot_scatter_from_df(
            df_train,
            out_png=figures_out / "Fig02E_acetate_train.png",
            title="Acetate-gradient (train anchors)",
            figsize=(6, 4),
        )
        rho_test = _plot_scatter_from_df(
            df_test,
            out_png=figures_out / "Fig02E_acetate_test.png",
            title="Acetate-gradient (test)",
            figsize=(6, 4),
        )
        (figures_out / "Fig02E_acetate_train_test_metrics.txt").write_text(
            f"rho_train={rho_train:.6g}\n" f"rho_test={rho_test:.6g}\n",
            encoding="utf-8",
        )
    except Exception as e:
        print(f"[WARN] Skipping Fig02E acetate train/test outputs: {e}")

    # Fig03/04: re-draw beeswarm with max_display=15 and figsize=(10,5)
    X_feat = _load_feature_matrix_regime_dataset(src_regime_dataset, seed=42, max_rows=2000)
    try:
        _plot_shap_beeswarm_from_saved_model(
            model_json=src_xgb_regime_model,
            X=X_feat,
            out_png=figures_out / "Fig03_shap_regime_beeswarm.png",
            is_classifier=True,
            max_display=15,
            figsize=(10, 5),
        )
        _plot_shap_beeswarm_from_saved_model(
            model_json=src_xgb_sev_model,
            X=X_feat,
            out_png=figures_out / "Fig04_shap_severity_beeswarm.png",
            is_classifier=False,
            max_display=15,
            figsize=(10, 5),
        )
    except RuntimeError as e:
        # Local dev environments may not have shap installed; keep script usable by falling back to re-saving PNGs.
        print(f"[WARN] Could not recompute SHAP beeswarm ({e}). Falling back to re-saving existing beeswarm PNGs.")
        _rewrite_png_with_figsize(src_shap_regime_png, figures_out / "Fig03_shap_regime_beeswarm.png", figsize=(10, 5))
        _rewrite_png_with_figsize(src_shap_sev_png, figures_out / "Fig04_shap_severity_beeswarm.png", figsize=(10, 5))

    # Fig05: copy DAG figure
    _copy(src_causal_png, figures_out / "Fig05_causal_dag.png")

    # Fig06: top edges bar (drawn)
    _plot_top_edges_bar(src_edge_stab, figures_out / "Fig06_top_edges_bar.png", topn=10)

    print(f"[OK] Wrote document: {doc_out}")
    print(f"[OK] Wrote figures under: {figures_out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

