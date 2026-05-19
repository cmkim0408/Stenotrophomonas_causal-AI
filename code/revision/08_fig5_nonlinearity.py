"""Workstream #D: Fig 5 non-linearity reanalysis.

Reviewer 2 questioned the "non-linear effects" claim around Fig 5b/5c. We
therefore re-derive SHAP dependence for the top-3 severity-regression
features (paper-named ADCS, APSR, APSR2) using OUT-OF-FOLD SHAP values from
5-fold CV (in-sample SHAP would leak) and apply a 3-test non-linearity
battery:

  (a) LOESS overlay (statsmodels.nonparametric.smoothers_lowess)
  (b) Piecewise-linear regression with grid-searched breakpoint
      (10–90 percentile), F-test piecewise vs single linear
  (c) Polynomial degree-2 OLS, F-test polynomial vs linear

For each feature we additionally compute SHAP interaction values
(Booster.predict pred_interactions=True) to surface top-K
interaction-modulated effects.

Verdict classification (verdict.txt):
  STRONG_NONLINEAR  — piecewise AND polynomial both p<0.05 AND |ΔR²|>0.10
  THRESHOLD_LIKE    — piecewise p<0.05 (breakpoint dominates) but
                      polynomial may or may not pass
  INTERACTION_MOD   — main-effect tests not significant but top SHAP
                      interaction value is sizeable
  LINEAR_RETREAT    — none of the tests pass

Outputs:
  revision_runs/iscience_rev1/08_nonlinearity/
    ├── fig5_dependence_with_fits.{png,pdf}
    ├── fig5_shap_interactions.{png,pdf}
    ├── nonlinearity_metrics.csv
    ├── nonlinearity_summary.md
    └── verdict.txt
"""
from __future__ import annotations

import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Sequence

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import xgboost as xgb
from scipy import stats
from sklearn.model_selection import KFold
from statsmodels.api import OLS, add_constant
from statsmodels.nonparametric.smoothers_lowess import lowess

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

OUT = u.ensure_outdirs("08_nonlinearity")
TOPK_FEATURES = 3
SEED = 42


@dataclass
class TestResult:
    feature: str
    test: str       # 'piecewise' | 'polynomial'
    p_value: float
    r2_baseline: float
    r2_model: float
    delta_r2: float
    extra: dict


# ---------------------------------------------------------------------------
# OOF SHAP via 5-fold CV
# ---------------------------------------------------------------------------

def compute_oof_shap(df: pd.DataFrame, feature_cols: list[str],
                     y: np.ndarray, *, n_splits: int = 5,
                     seed: int = SEED) -> tuple[np.ndarray, np.ndarray]:
    """Return (oof_preds, oof_shap) of shapes (n,), (n, n_feat)."""
    X = df[feature_cols].to_numpy()
    oof_pred = np.zeros(len(y), dtype=float)
    oof_shap = np.zeros((len(y), len(feature_cols)), dtype=float)
    cv = KFold(n_splits=n_splits, shuffle=True, random_state=seed)
    for fold, (tr, te) in enumerate(cv.split(X)):
        reg = xgb.XGBRegressor(
            n_estimators=300, max_depth=4, learning_rate=0.1,
            random_state=seed, n_jobs=-1, tree_method="hist", verbosity=0,
        )
        reg.fit(X[tr], y[tr])
        oof_pred[te] = reg.predict(X[te])
        contribs = reg.get_booster().predict(
            xgb.DMatrix(X[te]), pred_contribs=True)
        # contribs shape: (n_te, n_feat + 1); last col is bias
        oof_shap[te, :] = contribs[:, :-1]
    return oof_pred, oof_shap


# ---------------------------------------------------------------------------
# Non-linearity tests
# ---------------------------------------------------------------------------

def fit_linear(x: np.ndarray, y: np.ndarray) -> tuple[float, np.ndarray]:
    X = add_constant(x.reshape(-1, 1))
    res = OLS(y, X).fit()
    return float(res.rsquared), res.fittedvalues


def fit_polynomial(x: np.ndarray, y: np.ndarray, degree: int = 2
                   ) -> tuple[float, np.ndarray, float]:
    """Return (R², fitted, F-test p-value vs linear)."""
    n = len(x)
    X_lin = add_constant(x.reshape(-1, 1))
    X_pol = add_constant(np.column_stack([x, x**2]))
    res_lin = OLS(y, X_lin).fit()
    res_pol = OLS(y, X_pol).fit()
    # F-test: nested OLS comparison
    ssr_lin = float(res_lin.ssr)
    ssr_pol = float(res_pol.ssr)
    df_diff = res_lin.df_resid - res_pol.df_resid  # =1 for poly2 vs lin
    df_pol = float(res_pol.df_resid)
    if df_diff <= 0 or ssr_pol <= 0:
        return float(res_pol.rsquared), res_pol.fittedvalues, 1.0
    f_stat = ((ssr_lin - ssr_pol) / df_diff) / (ssr_pol / df_pol)
    p = float(1.0 - stats.f.cdf(f_stat, df_diff, df_pol))
    return float(res_pol.rsquared), res_pol.fittedvalues, p


def fit_piecewise_breakpoint(x: np.ndarray, y: np.ndarray
                             ) -> tuple[float, np.ndarray, float, float]:
    """Grid-search breakpoint between 10%-90% percentile; refit at best.
    Returns (best R², fitted at best, F-test p vs single linear, breakpoint).
    """
    lo, hi = np.percentile(x, 10), np.percentile(x, 90)
    candidates = np.linspace(lo, hi, 41)
    res_lin = OLS(y, add_constant(x.reshape(-1, 1))).fit()
    ssr_lin = float(res_lin.ssr)
    best = (0.0, None, 1.0, float(np.median(x)))
    for c in candidates:
        # Continuous piecewise: y = b0 + b1*x + b2*max(0, x-c)
        z = np.maximum(0.0, x - c)
        if z.var() == 0:  # all on one side
            continue
        X = add_constant(np.column_stack([x, z]))
        res = OLS(y, X).fit()
        ssr = float(res.ssr)
        df_diff = res_lin.df_resid - res.df_resid  # =1
        df_pw = float(res.df_resid)
        if df_diff <= 0 or ssr <= 0:
            continue
        f_stat = ((ssr_lin - ssr) / df_diff) / (ssr / df_pw)
        p = float(1.0 - stats.f.cdf(f_stat, df_diff, df_pw))
        r2 = float(res.rsquared)
        if r2 > best[0]:
            best = (r2, res.fittedvalues, p, float(c))
    if best[1] is None:
        best = (float(res_lin.rsquared), res_lin.fittedvalues, 1.0, float(np.median(x)))
    return best


# ---------------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------------

def plot_dependence_panel(ax, x: np.ndarray, y_shap: np.ndarray,
                          feature_short: str, piecewise_bp: float,
                          piecewise_pred: np.ndarray, poly_pred: np.ndarray,
                          loess_curve: np.ndarray, p_pw: float, p_poly: float):
    order = np.argsort(x)
    xs = x[order]
    # Scatter
    ax.scatter(x, y_shap, s=10, alpha=0.55, color="#0072B2",
               edgecolor="none", label="OOF SHAP")
    # LOESS
    ax.plot(loess_curve[:, 0], loess_curve[:, 1], color="#D55E00",
            lw=2.0, label="LOESS (frac=0.4)")
    # Piecewise
    ax.plot(xs, piecewise_pred[order], color="#009E73", lw=1.6,
            linestyle="--", label=f"piecewise (p={p_pw:.3g})")
    ax.axvline(piecewise_bp, color="#009E73", lw=0.8, linestyle=":",
               alpha=0.6, label=f"break = {piecewise_bp:.3g}")
    # Polynomial
    ax.plot(xs, poly_pred[order], color="#CC79A7", lw=1.4,
            linestyle="-.", label=f"poly2 (p={p_poly:.3g})")
    ax.set_xlabel(f"{feature_short} (FVA width)")
    ax.set_ylabel("OOF SHAP contribution to severity")
    ax.set_title(feature_short, fontsize=11)
    ax.grid(alpha=0.25)
    ax.legend(loc="best", fontsize=7.5, frameon=True)


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> None:
    print("[08 fig5 non-linearity] starting...")

    df = u.load_regime_dataset()
    feat_cols = u.width_cols(df)
    severity = u.severity_target(df).to_numpy()
    print(f"  dataset: {df.shape}  features: {len(feat_cols)}")

    # Step 1: get top-K severity SHAP features via OOF
    print("  Computing OOF SHAP on extended dataset (5-fold CV, seed=42)...")
    _oof_pred, oof_shap_full = compute_oof_shap(df, feat_cols, severity)
    feature_imp = np.abs(oof_shap_full).mean(axis=0)
    rank = np.argsort(-feature_imp)
    topk_idx = rank[:TOPK_FEATURES]
    top_features = [feat_cols[i] for i in topk_idx]
    print(f"  Top-{TOPK_FEATURES} severity features:")
    for i, fc in zip(topk_idx, top_features):
        print(f"    {fc:30s}  |SHAP|={feature_imp[i]:.5f}")

    # Step 2: per-feature 3-test battery on OOF SHAP
    test_rows: list[dict] = []
    fits = {}  # feature_col -> dict of curves for plotting
    for col_idx, col in zip(topk_idx, top_features):
        x = df[col].to_numpy().astype(float)
        y = oof_shap_full[:, col_idx].astype(float)
        # Linear baseline
        r2_lin, _ = fit_linear(x, y)
        # Polynomial
        r2_pol, poly_pred, p_poly = fit_polynomial(x, y, degree=2)
        # Piecewise
        r2_pw, pw_pred, p_pw, bp = fit_piecewise_breakpoint(x, y)
        # LOESS for plotting
        lo = lowess(y, x, frac=0.40, return_sorted=True)
        fits[col] = {
            "x": x, "y": y, "loess": lo,
            "piecewise_pred": pw_pred, "piecewise_bp": bp,
            "poly_pred": poly_pred,
            "p_pw": p_pw, "p_poly": p_poly,
            "r2_lin": r2_lin, "r2_pw": r2_pw, "r2_pol": r2_pol,
        }
        test_rows.extend([
            {"feature": col, "test": "linear_baseline",
             "p_value": float("nan"), "r2_baseline": r2_lin,
             "r2_model": r2_lin, "delta_r2": 0.0,
             "extra": ""},
            {"feature": col, "test": "piecewise_linear",
             "p_value": p_pw, "r2_baseline": r2_lin,
             "r2_model": r2_pw, "delta_r2": r2_pw - r2_lin,
             "extra": f"breakpoint={bp:.4f}"},
            {"feature": col, "test": "polynomial_d2",
             "p_value": p_poly, "r2_baseline": r2_lin,
             "r2_model": r2_pol, "delta_r2": r2_pol - r2_lin,
             "extra": ""},
        ])
        print(f"  {col:25s}  linR²={r2_lin:.3f}  pwR²={r2_pw:.3f} "
              f"(p={p_pw:.3g}, bp={bp:.3f})  polR²={r2_pol:.3f} (p={p_poly:.3g})")

    # Step 3: SHAP interactions (top-K pairwise)
    print("  Computing SHAP interaction values on a full-data fit...")
    reg = xgb.XGBRegressor(
        n_estimators=300, max_depth=4, learning_rate=0.1,
        random_state=SEED, n_jobs=-1, tree_method="hist", verbosity=0,
    )
    X = df[feat_cols].to_numpy()
    reg.fit(X, severity)
    interactions = reg.get_booster().predict(
        xgb.DMatrix(X), pred_interactions=True)
    # interactions shape: (n_samples, n_feat+1, n_feat+1)
    # Drop the bias row/col; symmetrize already
    inter = interactions[:, :-1, :-1]
    # Restrict to top-K features
    top_inter_panel = inter[:, topk_idx, :]  # (n, K, n_feat)
    # For each top feature, find best interacting partner (excluding self)
    interaction_partners = []
    for k, fi in enumerate(topk_idx):
        mean_abs_pair = np.abs(top_inter_panel[:, k, :]).mean(axis=0)
        mean_abs_pair[fi] = -1  # exclude self-interaction (main effect)
        partner_idx = int(np.argmax(mean_abs_pair))
        interaction_partners.append({
            "feature": top_features[k],
            "partner": feat_cols[partner_idx],
            "mean_abs_interaction": float(mean_abs_pair[partner_idx]),
            "partner_idx": partner_idx,
        })
    print("  Top SHAP interactions:")
    for p in interaction_partners:
        print(f"    {p['feature']:25s}  ←×→  {p['partner']:25s}  "
              f"|inter|={p['mean_abs_interaction']:.5f}")

    # ----- Verdict classification -----
    sig_pw = any(r["test"] == "piecewise_linear" and r["p_value"] < 0.05
                 for r in test_rows)
    sig_poly = any(r["test"] == "polynomial_d2" and r["p_value"] < 0.05
                   for r in test_rows)
    max_delta_pw = max((r["delta_r2"] for r in test_rows
                        if r["test"] == "piecewise_linear"), default=0.0)
    max_delta_pol = max((r["delta_r2"] for r in test_rows
                         if r["test"] == "polynomial_d2"), default=0.0)
    max_inter = max((p["mean_abs_interaction"] for p in interaction_partners),
                    default=0.0)
    # Heuristic thresholds per the prompt
    if sig_pw and sig_poly and max(max_delta_pw, max_delta_pol) > 0.10:
        verdict = "STRONG_NONLINEAR"
    elif sig_pw:
        verdict = "THRESHOLD_LIKE"
    elif max_inter > 0.005:   # ~ tenth of max main-effect SHAP
        verdict = "INTERACTION_MOD"
    else:
        verdict = "LINEAR_RETREAT"
    print(f"\n  >>> VERDICT: {verdict}")
    print(f"     piecewise sig: {sig_pw}; polynomial sig: {sig_poly}; "
          f"max Δ R² pw={max_delta_pw:.3f} poly={max_delta_pol:.3f}; "
          f"max |interaction|={max_inter:.5f}")

    # ----- Figure: dependence with fits (3 panels) -----
    fig, axes = plt.subplots(1, TOPK_FEATURES, figsize=(15, 4.4),
                             constrained_layout=True)
    if TOPK_FEATURES == 1:
        axes = [axes]
    for ax, col in zip(axes, top_features):
        fi = fits[col]
        short = col.replace("width__", "")
        plot_dependence_panel(ax, fi["x"], fi["y"], short,
                              fi["piecewise_bp"], fi["piecewise_pred"],
                              fi["poly_pred"], fi["loess"],
                              fi["p_pw"], fi["p_poly"])
    fig.suptitle(
        f"Fig 5 SHAP dependence with non-linearity diagnostics "
        f"(LOESS / piecewise / poly2 vs linear; OOF 5-fold; workstream #D). "
        f"Verdict = {verdict}.",
        fontsize=10.5)
    fig.savefig(OUT["ws"] / "fig5_dependence_with_fits.png", dpi=180,
                bbox_inches="tight")
    fig.savefig(OUT["ws"] / "fig5_dependence_with_fits.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote fig5_dependence_with_fits.{{png,pdf}}")

    # ----- Figure: SHAP interactions (top-3) -----
    fig, axes = plt.subplots(1, TOPK_FEATURES, figsize=(15, 4.4),
                             constrained_layout=True)
    if TOPK_FEATURES == 1:
        axes = [axes]
    for ax, p in zip(axes, interaction_partners):
        fi_col = p["feature"]
        partner_col = p["partner"]
        x = df[fi_col].to_numpy().astype(float)
        c = df[partner_col].to_numpy().astype(float)
        # Re-derive OOF SHAP for fi_col (use the cached oof_shap_full)
        col_idx_main = feat_cols.index(fi_col)
        y_shap = oof_shap_full[:, col_idx_main]
        # Colour by partner value
        sc = ax.scatter(x, y_shap, c=c, s=14, alpha=0.85,
                        cmap="viridis", edgecolor="none")
        cbar = fig.colorbar(sc, ax=ax, shrink=0.85)
        cbar.set_label(partner_col.replace("width__", ""),
                       fontsize=7.5)
        ax.set_xlabel(f"{fi_col.replace('width__','')} (FVA width)")
        ax.set_ylabel("OOF SHAP")
        ax.set_title(
            f"{fi_col.replace('width__','')}  modulated by  "
            f"{partner_col.replace('width__','')}\n"
            f"|interaction| = {p['mean_abs_interaction']:.4f}",
            fontsize=9.5)
        ax.grid(alpha=0.25)
    fig.suptitle(
        "SHAP interaction colouring: top-3 severity features "
        "vs their strongest interacting partners (workstream #D).",
        fontsize=10.5)
    fig.savefig(OUT["ws"] / "fig5_shap_interactions.png", dpi=180,
                bbox_inches="tight")
    fig.savefig(OUT["ws"] / "fig5_shap_interactions.pdf",
                bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote fig5_shap_interactions.{{png,pdf}}")

    # ----- nonlinearity_metrics.csv -----
    pd.DataFrame(test_rows).to_csv(
        OUT["ws"] / "nonlinearity_metrics.csv", index=False)
    pd.DataFrame(interaction_partners).to_csv(
        OUT["ws"] / "shap_interactions.csv", index=False)
    print("  wrote nonlinearity_metrics.csv + shap_interactions.csv")

    # ----- verdict.txt -----
    (OUT["ws"] / "verdict.txt").write_text(verdict + "\n", encoding="utf-8")
    print(f"  wrote verdict.txt: {verdict}")

    # ----- nonlinearity_summary.md -----
    md = ["# Fig 5 non-linearity reanalysis — summary",
          "",
          f"**Verdict:** `{verdict}`",
          "",
          "## What we did",
          "",
          ("Re-derived SHAP dependence for the top-3 severity-regression "
           f"features on the extended `regime_dataset_extended.parquet` "
           f"({len(df)} conditions × {len(feat_cols)} `width__` columns). "
           "Out-of-fold SHAP from 5-fold CV (random_state=42) was used to "
           "prevent leakage. Three non-linearity tests applied per feature:"),
          "",
          "- **LOESS overlay** — `statsmodels.nonparametric.smoothers_lowess` (frac=0.40)",
          "- **Piecewise-linear with grid-searched breakpoint** — best knot from 41-point "
          "grid between the 10th and 90th percentile; F-test piecewise vs single linear",
          "- **Polynomial degree-2** — OLS y ~ x + x² ; F-test vs linear",
          "",
          "## Top-3 severity features (mean |OOF SHAP|)",
          "",
          "| Rank | Feature | |OOF SHAP| |",
          "|---|---|---|"]
    for i, fc in zip(topk_idx, top_features):
        md.append(f"| {top_features.index(fc)+1} | `{fc}` | {feature_imp[i]:.5f} |")
    md += ["",
           "## Per-feature test results",
           "",
           "| Feature | Test | p-value | R² baseline (lin) | R² model | Δ R² | extra |",
           "|---|---|---|---|---|---|---|"]
    for r in test_rows:
        md.append(
            f"| `{r['feature']}` | {r['test']} | {r['p_value']:.4g} | "
            f"{r['r2_baseline']:.3f} | {r['r2_model']:.3f} | "
            f"{r['delta_r2']:+.3f} | {r['extra']} |")
    md += ["",
           "## SHAP interaction partners (top-3 features)",
           "",
           "| Feature | Strongest interacting partner | mean &#124;interaction&#124; |",
           "|---|---|---|"]
    for p in interaction_partners:
        md.append(
            f"| `{p['feature']}` | `{p['partner']}` | {p['mean_abs_interaction']:.5f} |")
    md += ["",
           "## Verdict logic",
           "",
           f"- Any piecewise test with p<0.05? **{sig_pw}**",
           f"- Any polynomial test with p<0.05? **{sig_poly}**",
           f"- Max Δ R² over linear: piecewise = {max_delta_pw:+.3f}, "
           f"polynomial = {max_delta_pol:+.3f}",
           f"- Max mean |interaction| (top-3 features) = {max_inter:.5f}",
           "",
           f"→ Verdict = `{verdict}`",
           "",
           "Branching rule (per the workstream-D prompt):",
           "",
           "  - `STRONG_NONLINEAR` if (piecewise sig.) AND (polynomial sig.) AND (max Δ R² > 0.10)",
           "  - `THRESHOLD_LIKE`   if piecewise sig. (independent of polynomial)",
           "  - `INTERACTION_MOD`  if main-effect tests not sig. but max |interaction| > 0.005",
           "  - `LINEAR_RETREAT`   otherwise",
           "",
           "## Recommended caption variant",
           "",
           "See `docs/revision/captions/fig5_caption_v2.md`. The variant marked "
           f"`RECOMMENDED based on workstream #D verdict: {verdict}` should be "
           "used for the Fig 5 caption (or the SI dependence-plot caption if "
           "main Fig 5 is left untouched).",
           "",
           "## Honest interpretation",
           ""]
    if verdict == "STRONG_NONLINEAR":
        md.append("Both the piecewise-linear and polynomial fits significantly "
                  "improve over the linear baseline (p<0.05) with Δ R² > 0.10. "
                  "The 'non-linear effects' wording in the manuscript is "
                  "supported. The breakpoint provides a data-driven threshold "
                  "for the flexibility-collapse interpretation.")
    elif verdict == "THRESHOLD_LIKE":
        md.append("The piecewise-linear fit significantly improves over the "
                  "linear baseline (p<0.05) with a clear breakpoint, supporting "
                  "a 'threshold-like' interpretation. The polynomial test may "
                  "or may not reach significance, which is consistent with a "
                  "discontinuous change in slope rather than a smooth curve.")
    elif verdict == "INTERACTION_MOD":
        md.append("Main-effect non-linearity tests do not reach significance, "
                  "but a sizeable SHAP-interaction value indicates that the "
                  "dependence is modulated by a partner reaction. The "
                  "'non-linear effects' wording in the manuscript is more "
                  "accurately rephrased as 'interaction-modulated effects'.")
    else:  # LINEAR_RETREAT
        md.append("None of the three non-linearity tests reach significance "
                  "within the tested LHS range. The 'non-linear effects' "
                  "wording in the manuscript should be revised to the more "
                  "conservative 'monotonic trends consistent with progressive "
                  "rigidification within the tested condition envelope'. "
                  "Stronger non-linearity may emerge at constraint regimes "
                  "outside the sampled LHS.")
    (OUT["ws"] / "nonlinearity_summary.md").write_text(
        "\n".join(md) + "\n", encoding="utf-8")
    print(f"  wrote nonlinearity_summary.md")

    print("[08 fig5 non-linearity] done.")


if __name__ == "__main__":
    main()
