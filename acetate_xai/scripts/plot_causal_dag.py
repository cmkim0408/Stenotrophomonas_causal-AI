from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd


EXOGENOUS = {"atpm_fixed", "o2_lb", "acetate_mode", "nh4_mode", "frac_opt", "targets_n"}
OUTCOMES = {"primary_regime", "maintenance_severity"}


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Plot causal DAG from edge list + stability table.")
    p.add_argument("--edges", required=True, help="dag_edges.csv path")
    p.add_argument("--stability", required=True, help="edge_stability.csv path")
    p.add_argument("--outdir", required=True, help="Output directory (e.g., results/causal_dag)")
    p.add_argument("--seed", type=int, default=42, help="Layout random seed")
    p.add_argument("--verbose", "-v", action="store_true", help="Verbose logging")
    return p


def _edge_key(src: str, dst: str, kind: str) -> str:
    if kind == "directed":
        return f"{src}->{dst}"
    a, b = sorted([src, dst])
    return f"{a}--{b}"


def _safe_imports():
    try:
        import networkx as nx  # noqa: F401
    except Exception as e:  # noqa: BLE001
        raise RuntimeError(f"Missing dependency: networkx ({e}). Install: pip install networkx") from e

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt  # noqa: F401


def _node_color(name: str) -> str:
    if name in EXOGENOUS:
        return "#4C78A8"  # blue
    if name in OUTCOMES:
        return "#E45756"  # red
    return "#72B7B2"  # teal (SHAP features / others)


def _sanitize_label(s: str) -> str:
    # For nicer plots: shorten very long feature names
    if len(s) <= 35:
        return s
    return s[:32] + "..."


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    _safe_imports()

    import networkx as nx
    import matplotlib.pyplot as plt

    edges_path = Path(args.edges)
    stab_path = Path(args.stability)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    if not edges_path.exists():
        raise FileNotFoundError(f"--edges not found: {edges_path}")
    if not stab_path.exists():
        raise FileNotFoundError(f"--stability not found: {stab_path}")

    edges = pd.read_csv(edges_path)
    if not {"src", "dst", "kind"}.issubset(edges.columns):
        raise ValueError("dag_edges.csv must have columns: src, dst, kind")

    stab = pd.read_csv(stab_path)
    if "edge" not in stab.columns or "frequency" not in stab.columns:
        raise ValueError("edge_stability.csv must have columns: edge, frequency")
    stab_map = dict(zip(stab["edge"].astype(str), pd.to_numeric(stab["frequency"], errors="coerce").fillna(0.0)))

    # Build graph
    G = nx.DiGraph()
    for _, r in edges.iterrows():
        src = str(r["src"])
        dst = str(r["dst"])
        kind = str(r["kind"])
        k = _edge_key(src, dst, kind)
        freq = float(stab_map.get(k, 0.0))
        G.add_node(src)
        G.add_node(dst)
        if kind == "directed":
            G.add_edge(src, dst, kind=kind, frequency=freq)
        else:
            # store as undirected-like edge, but draw as two-way dashed in DiGraph
            G.add_edge(src, dst, kind=kind, frequency=freq)

    # Layout
    pos = nx.spring_layout(G, seed=int(args.seed))

    # Node styling
    nodes = list(G.nodes())
    node_colors = [_node_color(n) for n in nodes]
    node_sizes = [900 if (n in EXOGENOUS or n in OUTCOMES) else 600 for n in nodes]

    # Edge styling based on stability frequency
    edges_list = list(G.edges(data=True))
    widths = [1.0 + 4.0 * float(d.get("frequency", 0.0)) for _, _, d in edges_list]
    alphas = [0.2 + 0.8 * float(d.get("frequency", 0.0)) for _, _, d in edges_list]

    fig, ax = plt.subplots(figsize=(12, 8))
    ax.set_axis_off()

    nx.draw_networkx_nodes(G, pos, node_color=node_colors, node_size=node_sizes, ax=ax, linewidths=1.0, edgecolors="#333333")

    # Draw directed edges with per-edge alpha by splitting into segments
    for (u, v, d), w, a in zip(edges_list, widths, alphas):
        style = "-" if d.get("kind") == "directed" else "--"
        nx.draw_networkx_edges(
            G,
            pos,
            edgelist=[(u, v)],
            width=w,
            alpha=a,
            arrows=True,
            arrowstyle="-|>",
            arrowsize=16,
            edge_color="#333333",
            style=style,
            ax=ax,
        )

    labels = {n: _sanitize_label(n) for n in nodes}
    nx.draw_networkx_labels(G, pos, labels=labels, font_size=8, ax=ax)

    # Legend (proxy artists)
    from matplotlib.lines import Line2D

    legend_elems = [
        Line2D([0], [0], marker="o", color="w", label="Exogenous", markerfacecolor="#4C78A8", markersize=10),
        Line2D([0], [0], marker="o", color="w", label="SHAP feature", markerfacecolor="#72B7B2", markersize=10),
        Line2D([0], [0], marker="o", color="w", label="Outcome", markerfacecolor="#E45756", markersize=10),
    ]
    ax.legend(handles=legend_elems, loc="lower left", fontsize=9, frameon=True)

    out_png = outdir / "causal_dag.png"
    fig.savefig(out_png, dpi=220, bbox_inches="tight")
    print(f"[OK] Wrote: {out_png}")

    # PDF (best-effort)
    out_pdf = outdir / "causal_dag.pdf"
    try:
        fig.savefig(out_pdf, bbox_inches="tight")
        print(f"[OK] Wrote: {out_pdf}")
    except Exception:
        pass
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

