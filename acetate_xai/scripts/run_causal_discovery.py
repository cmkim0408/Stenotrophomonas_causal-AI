from __future__ import annotations

import argparse
import json
import logging
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd


@dataclass(frozen=True)
class Edge:
    src: str
    dst: str
    kind: str  # "directed" | "undirected"

    def key(self) -> str:
        if self.kind == "directed":
            return f"{self.src}->{self.dst}"
        return "--".join(sorted([self.src, self.dst]))


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Run causal discovery (PC + bootstrap) on causal_dataset.")
    p.add_argument("--data", default="results/causal_dataset.parquet", help="Input causal dataset parquet")
    p.add_argument("--topk-regime", required=True, help="CSV with regime top features (expects column 'feature')")
    p.add_argument("--topk-severity", required=True, help="CSV with severity top features (expects column 'feature')")
    p.add_argument("--outdir", default="results/causal_dag", help="Output directory")
    p.add_argument("--bootstrap", type=int, default=100, help="Bootstrap iterations")
    p.add_argument("--seed", type=int, default=42, help="Random seed")
    p.add_argument("--method", default="pc", choices=["pc", "notears"], help="Discovery method")
    p.add_argument("--use-priors", action="store_true", help="Apply BackgroundKnowledge priors (default: off).")
    p.add_argument(
        "--indep-test",
        default="fisherz",
        choices=["fisherz", "kci"],
        help="PC independence test (default: fisherz). If fisherz fails, kci will be attempted as fallback when available.",
    )
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _feature_cols(df: pd.DataFrame) -> list[str]:
    return [c for c in df.columns if c.startswith(("width__", "mid__", "signchange__"))]


def _load_topk(path: str | Path, k: int) -> list[str]:
    df = pd.read_csv(path)
    if "feature" not in df.columns:
        raise ValueError(f"topk file missing column 'feature': {path}")
    feats = df["feature"].astype(str).tolist()
    return [f for f in feats if f][:k]


def _encode_categoricals(df: pd.DataFrame, cols: Iterable[str]) -> pd.DataFrame:
    out = df.copy()
    for c in cols:
        if c in out.columns:
            out[c] = out[c].astype("category").cat.codes.replace({-1: np.nan})
    return out


def _remove_constant_columns(df: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    """
    Remove columns with variance == 0 (including all-constant columns).
    """
    # ddof=0 gives population variance; robust for small samples
    var = df.var(axis=0, ddof=0)
    const_cols = var.index[(var == 0) | var.isna()].astype(str).tolist()
    if const_cols:
        return df.drop(columns=const_cols), const_cols
    return df, []


def _remove_collinear_columns(df: pd.DataFrame, *, threshold: float = 0.9999) -> tuple[pd.DataFrame, list[str]]:
    """
    Remove columns to avoid near-perfect collinearity that makes correlation matrices singular.

    Strategy: compute abs(corr), look at upper triangle, and drop the second column in any pair
    with abs(corr) >= threshold.
    """
    if df.shape[1] <= 1:
        return df, []

    corr = df.corr().abs()
    # upper triangle mask (excluding diagonal)
    upper = corr.where(np.triu(np.ones(corr.shape), k=1).astype(bool))
    to_drop: list[str] = []
    for col in upper.columns:
        if (upper[col] >= threshold).any():
            to_drop.append(str(col))
    if to_drop:
        return df.drop(columns=to_drop), to_drop
    return df, []


def _add_jitter(df: pd.DataFrame, *, seed: int, scale: float) -> pd.DataFrame:
    """
    Add tiny noise to avoid singular correlation matrices in FisherZ:
      noise = scale * std(col) * N(0,1)
    """
    if df.empty:
        return df
    rng = np.random.default_rng(int(seed))
    std = df.std(axis=0, ddof=0).replace({0.0: np.nan}).fillna(1.0).to_numpy()
    noise = rng.normal(0.0, 1.0, size=df.shape) * (float(scale) * std)
    out = df.to_numpy(dtype=float) + noise
    return pd.DataFrame(out, columns=df.columns, index=df.index)


def _pc_preprocess(
    work: pd.DataFrame,
    *,
    seed: int,
    corr_threshold: float = 0.9999,
) -> tuple[pd.DataFrame, dict]:
    """
    Preprocess before PC(FisherZ) to prevent singular correlation matrix:
    - numeric only coercion
    - remove variance==0 cols
    - remove collinear cols (abs corr >= threshold)
    """
    df = work.copy()
    # force numeric
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.fillna(0.0)

    before = df.shape[1]
    df, removed_constant = _remove_constant_columns(df)
    df, removed_collinear = _remove_collinear_columns(df, threshold=corr_threshold)
    meta = {
        "removed_constant": removed_constant,
        "removed_collinear": removed_collinear,
        "final_vars": int(df.shape[1]),
        "initial_vars": int(before),
        "corr_threshold": float(corr_threshold),
        "jitter_base": 1e-8,
        "jitter_retry": 1e-7,
        "seed": int(seed),
    }
    return df, meta


def _prepare_causal_matrix(
    df: pd.DataFrame,
    variables: list[str],
    *,
    min_rows: int = 200,
) -> tuple[pd.DataFrame, dict]:
    """
    Prepare a numeric matrix for causal discovery without dropping almost all rows.

    Rules (minimal change, robust):
    - SHAP feature columns: numeric, fillna(0.0)
    - Numeric meta columns (atpm_fixed, o2_lb, frac_opt, targets_n): create has_<col>, fillna(0.0)
    - Categorical meta columns (acetate_mode, nh4_mode): create has_<col>, fillna("none"), one-hot encode
    - Outcomes: keep; encode primary_regime to codes; drop rows only if outcomes are missing
    """
    work = df[variables].copy()

    # Identify groups
    shap_feats = [c for c in work.columns if c.startswith(("width__", "mid__", "signchange__"))]
    num_meta = [c for c in ["atpm_fixed", "o2_lb", "frac_opt", "targets_n"] if c in work.columns]
    cat_meta = [c for c in ["acetate_mode", "nh4_mode"] if c in work.columns]
    outcomes = [c for c in ["primary_regime", "maintenance_severity"] if c in work.columns]

    meta_info = {"shap_feats": shap_feats, "num_meta": num_meta, "cat_meta": cat_meta, "outcomes": outcomes}

    # Outcomes: we will drop rows only if outcomes are missing
    if "primary_regime" in work.columns:
        work["primary_regime"] = work["primary_regime"].astype(str)
        # restore NaN-like strings to NaN
        work.loc[work["primary_regime"].str.strip().str.lower().isin(["", "nan"]), "primary_regime"] = np.nan
        work["primary_regime"] = work["primary_regime"].astype("category").cat.codes.replace({-1: np.nan})
    if "maintenance_severity" in work.columns:
        work["maintenance_severity"] = pd.to_numeric(work["maintenance_severity"], errors="coerce")

    # Numeric meta: has_ indicator + fillna(0)
    for c in num_meta:
        work[f"has_{c}"] = (~pd.to_numeric(work[c], errors="coerce").isna()).astype(int)
        work[c] = pd.to_numeric(work[c], errors="coerce").fillna(0.0)

    # Categorical meta: has_ indicator + fillna("none") + one-hot
    for c in cat_meta:
        has = work[c].notna() & (work[c].astype(str).str.strip() != "") & (work[c].astype(str).str.lower() != "nan")
        work[f"has_{c}"] = has.astype(int)
        work[c] = work[c].astype(str)
        work.loc[~has, c] = "none"
        dummies = pd.get_dummies(work[c], prefix=c, prefix_sep="=", dtype=int)
        work = pd.concat([work.drop(columns=[c]), dummies], axis=1)

    # SHAP features: numeric fill
    for c in shap_feats:
        work[c] = pd.to_numeric(work[c], errors="coerce").fillna(0.0)

    # Final: drop rows only if outcomes are missing (if present)
    if outcomes:
        work = work.dropna(subset=outcomes).copy()

    # Ensure numeric everywhere
    for c in work.columns:
        work[c] = pd.to_numeric(work[c], errors="coerce")
    work = work.fillna(0.0)

    if len(work) < min_rows:
        raise RuntimeError(
            f"Too few rows after meta-aware cleaning: {len(work)} (need >= {min_rows}). "
            "Check whether outcomes are missing or the dataset is too small."
        )

    return work, meta_info


def _select_variables(
    df: pd.DataFrame,
    *,
    top_regime_15: list[str],
    top_sev_15: list[str],
) -> tuple[list[str], dict]:
    """
    Build variable list:
      exogenous (if exists) + shap topk (dedup) + outcomes
    With fallback if variable count > 30: reduce to top10 + top10 (keep exogenous+outcomes).
    """
    exogenous_candidates = ["atpm_fixed", "o2_lb", "acetate_mode", "nh4_mode", "frac_opt", "targets_n"]
    exogenous = [c for c in exogenous_candidates if c in df.columns]

    # SHAP top-k variables: must exist in df and be true feature columns
    feat_set = set(_feature_cols(df))
    shap_union_15 = []
    for f in (top_regime_15 + top_sev_15):
        if f in feat_set and f not in shap_union_15:
            shap_union_15.append(f)

    outcomes = []
    if "primary_regime" in df.columns:
        outcomes.append("primary_regime")
    if "maintenance_severity" in df.columns:
        outcomes.append("maintenance_severity")

    vars_all = exogenous + shap_union_15 + outcomes

    # fallback
    if len(vars_all) > 30:
        top_regime_10 = [f for f in top_regime_15[:10] if f in feat_set]
        top_sev_10 = [f for f in top_sev_15[:10] if f in feat_set]
        shap_union_10 = []
        for f in (top_regime_10 + top_sev_10):
            if f not in shap_union_10:
                shap_union_10.append(f)
        vars_all = exogenous + shap_union_10 + outcomes

    meta = {
        "n_exogenous": len(exogenous),
        "n_shap_features": len(vars_all) - len(exogenous) - len(outcomes),
        "n_outcomes": len(outcomes),
        "n_total": len(vars_all),
        "exogenous": exogenous,
        "outcomes": outcomes,
    }
    return vars_all, meta


def _apply_priors_pc(background_knowledge, variables: list[str]) -> dict:
    """
    Priors:
    - exogenous variables: forbid any incoming edges into them
    - maintenance_severity: forbid outgoing edges from it
    """
    priors = {"forbid_incoming_to_exogenous": [], "forbid_outgoing_from_maintenance": []}
    # After preprocessing, exogenous variables can expand to has_* and one-hot columns.
    exo_prefixes = ("atpm_fixed", "o2_lb", "frac_opt", "targets_n", "acetate_mode", "nh4_mode", "has_")
    exogenous = [v for v in variables if v.startswith(exo_prefixes)]
    maint = "maintenance_severity" if "maintenance_severity" in variables else None

    # forbid X -> exog for all X != exog
    for ex in exogenous:
        for src in variables:
            if src == ex:
                continue
            background_knowledge.add_forbidden_by_node(src, ex)
            priors["forbid_incoming_to_exogenous"].append(f"{src}->{ex}")

    # forbid maintenance -> others
    if maint:
        for dst in variables:
            if dst == maint:
                continue
            background_knowledge.add_forbidden_by_node(maint, dst)
            priors["forbid_outgoing_from_maintenance"].append(f"{maint}->{dst}")

    return priors


def _edges_from_causallearn_graph(g, variables: list[str]) -> list[Edge]:
    """
    Convert causal-learn graph to edges. We emit directed edges when orientation is present,
    otherwise undirected (skeleton).
    """
    # causal-learn uses endpoints codes; easiest is to inspect graph endpoints
    # g.G.graph is a matrix of endpoint codes (N x N)
    mat = np.asarray(g.G.graph)
    n = len(variables)
    edges: list[Edge] = []
    for i in range(n):
        for j in range(i + 1, n):
            a = mat[i, j]
            b = mat[j, i]
            # Endpoint codes per causal-learn:
            # 0: NULL, 1: ARROW, 2: TAIL, 3: CIRCLE (varies). We'll handle common cases.
            # Directed i -> j typically: mat[i,j]=ARROW and mat[j,i]=TAIL
            if (a == 1 and b == 2) or (a == 1 and b == 3):  # arrow at j, tail/circle at i
                edges.append(Edge(src=variables[i], dst=variables[j], kind="directed"))
            elif (b == 1 and a == 2) or (b == 1 and a == 3):
                edges.append(Edge(src=variables[j], dst=variables[i], kind="directed"))
            elif a != 0 or b != 0:
                # some connection exists but not confidently directed
                edges.append(Edge(src=variables[i], dst=variables[j], kind="undirected"))
    return edges


def _run_pc(
    data: np.ndarray,
    variables: list[str],
    seed: int,
    use_priors: bool,
    indep_test: str,
) -> tuple[list[Edge], dict]:
    """
    PC Algorithm using causal-learn (recommended).
    """
    try:
        from causallearn.search.ConstraintBased.PC import pc
        from causallearn.utils.PCUtils.BackgroundKnowledge import BackgroundKnowledge
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            "Missing dependency for PC: causal-learn. Install on server:\n"
            "  pip install causal-learn\n"
            f"Import error: {e}"
        ) from e

    bk = None
    priors: dict = {}
    if use_priors:
        bk = BackgroundKnowledge()
        priors = _apply_priors_pc(bk, variables)
        logging.info("[INFO] priors enabled (--use-priors)")
    else:
        logging.info("[INFO] priors skipped (default)")

    # Use FisherZ by default. If it fails due to singular correlation matrices,
    # caller is expected to jitter/retry; we also allow kci when requested/available.
    g = pc(data, alpha=0.05, indep_test=indep_test, stable=True, uc_rule=0, uc_priority=2, background_knowledge=bk)
    edges = _edges_from_causallearn_graph(g, variables)
    meta = {"alpha": 0.05, "indep_test": indep_test, "priors": priors, "seed": seed, "use_priors": use_priors}
    return edges, meta


def _run_notears(data: np.ndarray, variables: list[str], seed: int) -> tuple[list[Edge], dict]:
    """
    Optional NOTEARS (if dependency available). Implemented as a best-effort fallback.
    """
    try:
        from notears import linear  # type: ignore
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(
            "NOTEARS selected but dependency is missing. Install on server:\n"
            "  pip install notears\n"
            f"Import error: {e}"
        ) from e

    # NOTEARS expects continuous; we use the numeric matrix directly.
    w_est = linear.notears_linear(data, lambda1=0.1, loss_type="l2", max_iter=100, w_threshold=0.3)
    edges: list[Edge] = []
    n = len(variables)
    for i in range(n):
        for j in range(n):
            if i == j:
                continue
            if abs(w_est[i, j]) > 0:
                edges.append(Edge(src=variables[i], dst=variables[j], kind="directed"))
    meta = {"seed": seed, "lambda1": 0.1, "w_threshold": 0.3}
    return edges, meta


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(message)s")

    data_path = Path(args.data)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        print(f"[ERROR] data not found: {data_path}")
        return 2

    df = pd.read_parquet(data_path)

    # Load SHAP top-k feature lists
    top_regime_15 = _load_topk(args.topk_regime, k=15)
    top_sev_15 = _load_topk(args.topk_severity, k=15)

    variables, sel_meta = _select_variables(df, top_regime_15=top_regime_15, top_sev_15=top_sev_15)
    if len(variables) < 3:
        print("[ERROR] Not enough variables selected for causal discovery (need >= 3).")
        return 2

    # Prepare numeric data matrix (meta-aware; do NOT drop rows just because meta is NaN)
    try:
        work, prep_meta = _prepare_causal_matrix(df, variables, min_rows=200)
    except Exception as e:  # noqa: BLE001
        print(f"[ERROR] {e}")
        return 2

    rng = np.random.default_rng(int(args.seed))

    # Preprocess for PC to avoid singular correlation matrices
    work_pc, pc_meta = _pc_preprocess(work, seed=int(args.seed))
    print(f"[INFO] removed_constant={len(pc_meta['removed_constant'])} removed_collinear={len(pc_meta['removed_collinear'])} final_vars={pc_meta['final_vars']}")

    # Fit once on full cleaned data
    if args.method == "pc":
        # Jitter strategy for FisherZ singular issue
        indep = str(args.indep_test)
        try:
            data_mat = _add_jitter(work_pc, seed=int(args.seed), scale=pc_meta["jitter_base"]).to_numpy(dtype=float)
            edges, method_meta = _run_pc(
                data_mat,
                variables=list(work_pc.columns),
                seed=int(args.seed),
                use_priors=bool(args.use_priors),
                indep_test=indep,
            )
        except Exception as e1:  # noqa: BLE001
            logging.warning("PC failed on first attempt (%s). Retrying with stronger jitter...", e1)
            try:
                data_mat = _add_jitter(work_pc, seed=int(args.seed), scale=pc_meta["jitter_retry"]).to_numpy(dtype=float)
                edges, method_meta = _run_pc(
                    data_mat,
                    variables=list(work_pc.columns),
                    seed=int(args.seed),
                    use_priors=bool(args.use_priors),
                    indep_test=indep,
                )
            except Exception as e2:  # noqa: BLE001
                if indep != "kci":
                    logging.warning("PC still failing with fisherz. Trying indep_test=kci as fallback...")
                    data_mat = work_pc.to_numpy(dtype=float)
                    edges, method_meta = _run_pc(
                        data_mat,
                        variables=list(work_pc.columns),
                        seed=int(args.seed),
                        use_priors=bool(args.use_priors),
                        indep_test="kci",
                    )
                else:
                    raise
        method_meta["pc_preprocess"] = pc_meta
    else:
        data_mat = work_pc.to_numpy(dtype=float)
        edges, method_meta = _run_notears(data_mat, variables=list(work_pc.columns), seed=int(args.seed))

    # Bootstrap stability
    B = int(args.bootstrap)
    if B <= 0:
        print("[ERROR] --bootstrap must be > 0")
        return 2

    edge_counts: dict[str, int] = {}
    for b in range(B):
        idx = rng.integers(0, len(work), size=len(work))
        boot_df = work.iloc[idx].copy()
        boot_df, _ = _pc_preprocess(boot_df, seed=int(args.seed) + b)
        boot = _add_jitter(boot_df, seed=int(args.seed) + b, scale=pc_meta["jitter_base"]).to_numpy(dtype=float)
        try:
            if args.method == "pc":
                # same indep_test setting as main run
                e_b, _ = _run_pc(
                    boot,
                    variables=list(boot_df.columns),
                    seed=int(args.seed),
                    use_priors=bool(args.use_priors),
                    indep_test=str(args.indep_test),
                )
            else:
                e_b, _ = _run_notears(boot, variables=list(boot_df.columns), seed=int(args.seed))
        except Exception as e:  # noqa: BLE001
            logging.warning("Bootstrap iteration %d failed: %s", b, e)
            continue
        keys = {e.key() for e in e_b}
        for k in keys:
            edge_counts[k] = edge_counts.get(k, 0) + 1

    # Outputs
    dag_edges = pd.DataFrame([asdict(e) for e in edges])
    dag_edges.to_csv(outdir / "dag_edges.csv", index=False)

    stab = pd.DataFrame(
        [{"edge": k, "frequency": v / B, "count": v, "bootstrap": B} for k, v in sorted(edge_counts.items(), key=lambda x: -x[1])]
    )
    stab.to_csv(outdir / "edge_stability.csv", index=False)

    meta = {
        "data": str(data_path),
        "method": args.method,
        "seed": int(args.seed),
        "bootstrap": B,
        "variables": list(work.columns),
        "selection": sel_meta,
        "method_meta": method_meta,
        "priors": method_meta.get("priors", {}),
        "n_rows_used": int(len(work)),
    }
    (outdir / "dag_metadata.json").write_text(json.dumps(meta, indent=2), encoding="utf-8")

    print(f"[INFO] variables={len(work.columns)} rows_used={len(work)} method={args.method} bootstrap={B}")
    print(f"[OK] Wrote: {outdir / 'dag_edges.csv'} (edges={len(dag_edges)})")
    print(f"[OK] Wrote: {outdir / 'edge_stability.csv'} (unique_edges={len(stab)})")
    print(f"[OK] Wrote: {outdir / 'dag_metadata.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

