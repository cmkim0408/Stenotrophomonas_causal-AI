"""Workstream #E: Fig 6 arrow consistency audit.

The PC bootstrap (acetate_xai/scripts/run_causal_discovery.py) emits
*undirected* edges into results/figures_final/data/Fig05_causal_dag__dag_edges.csv
(all rows have kind = 'undirected'); the corresponding bootstrap stability scores
in Fig05_causal_dag__edge_stability.csv use the '--' separator (also undirected).

The published Fig 6a (results/figures_final/Fig05_causal_dag.png) nonetheless
renders arrowheads, and Fig 6b labels use the '->' arrow notation. Reviewer 2
flagged this representation inconsistency.

The honest fix is to re-render both panels with a consistent UNDIRECTED
convention that matches the underlying algorithmic output, and to explicitly
state in the caption that the PC algorithm produced an equivalence class of
undirected edges (no orientation rule succeeded within the bootstrap-stable set).

Outputs:
  revision_runs/iscience_rev1/09_fig6_audit/fig6_consistent.{png,pdf,svg}
  revision_runs/iscience_rev1/09_fig6_audit/fig6_before_after.{png,pdf}
  revision_runs/iscience_rev1/09_fig6_audit/fig6_audit_summary.md
"""
from __future__ import annotations

import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

OUT_BASE = u.REVISION_OUT / "09_fig6_audit"
OUT_BASE.mkdir(parents=True, exist_ok=True)

EDGE_STAB_CSV = u.ROOT / "results" / "figures_final" / "data" / "Fig05_causal_dag__edge_stability.csv"
DAG_EDGES_CSV = u.ROOT / "results" / "figures_final" / "data" / "Fig05_causal_dag__dag_edges.csv"
PUBLISHED_6A = u.ROOT / "results" / "figures_final" / "Fig05_causal_dag.png"
PUBLISHED_6B = u.ROOT / "results" / "figures_draft" / "Fig05b_causal_stability_clean.png"

# Module-level readable label map (shortened paper-named conventions)
READABLE = {
    "atpm_fixed": "ATPM (fixed)",
    "o2_lb": "O$_2$ uptake bound",
    "primary_regime": "Regime (O$_2$/N/Ac)",
    "maintenance_severity": "Growth severity",
    "width__ACONT": "Aconitase Flex.",
    "width__ADCS": "ADC Synthase Flex.",
    "width__AKGDH": "AKGDH Flex.",
    "width__CYO1_KT": "CYO1_KT Flex.",
    "width__DMPPS": "DMPPS Flex.",
    "width__EX_h_e": "EX_H$^+$ Flex.",
    "width__EX_h2o_e": "EX_H$_2$O Flex.",
    "width__EX_o2_e": "EX_O$_2$ Flex.",
    "width__EX_co2_e": "EX_CO$_2$ Flex.",
    "width__12DGR120tipp": "1,2-DAG transport Flex.",
    "width__5DOAN": "5-DOA synthase Flex.",
    "mid__EX_h2o_e": "EX_H$_2$O Mid.",
    "mid__EX_o2_e": "EX_O$_2$ Mid.",
    "mid__EX_h_e": "EX_H$^+$ Mid.",
    "mid__AHCi": "AHCi Mid.",
}


def short(node: str) -> str:
    return READABLE.get(node, node)


def _parse_edge(edge_str: str) -> tuple[str, str]:
    if "--" in edge_str:
        a, b = edge_str.split("--", 1)
        return a.strip(), b.strip()
    if "->" in edge_str:
        a, b = edge_str.split("->", 1)
        return a.strip(), b.strip()
    raise ValueError(f"Unrecognized edge: {edge_str!r}")


def load_data() -> pd.DataFrame:
    stab = pd.read_csv(EDGE_STAB_CSV)
    parsed = stab["edge"].apply(_parse_edge)
    stab["src"] = parsed.apply(lambda t: t[0])
    stab["dst"] = parsed.apply(lambda t: t[1])
    return stab


def draw_undirected_graph(ax, edges: pd.DataFrame, top_k: int = 12,
                          seed: int = 42, title: str = "(a) Undirected adjacency from PC bootstrap"):
    edges = edges.sort_values("frequency", ascending=False).head(top_k).copy()
    G = nx.Graph()
    for _, r in edges.iterrows():
        G.add_edge(r["src"], r["dst"], frequency=float(r["frequency"]))
    pos = nx.spring_layout(G, k=1.6, seed=seed, iterations=200)

    # Edge widths proportional to bootstrap frequency
    freqs = np.array([G[u][v]["frequency"] for u, v in G.edges()])
    widths = 0.8 + 4.0 * freqs

    # Node colours: outcomes vs exogenous vs feature
    node_colors = []
    for n in G.nodes():
        if n in ("primary_regime", "maintenance_severity"):
            node_colors.append("#D55E00")
        elif n in ("o2_lb", "atpm_fixed"):
            node_colors.append("#0072B2")
        else:
            node_colors.append("#56B4E9")

    nx.draw_networkx_edges(G, pos, ax=ax, width=widths,
                           edge_color="#555555", alpha=0.75,
                           arrows=False)  # ← KEY: undirected, no arrowheads
    nx.draw_networkx_nodes(G, pos, ax=ax, node_color=node_colors,
                           node_size=750, edgecolors="black", linewidths=0.7)
    nx.draw_networkx_labels(G, pos, ax=ax,
                            labels={n: short(n) for n in G.nodes()},
                            font_size=7.5)

    # Legend dummies
    from matplotlib.lines import Line2D
    handles = [
        Line2D([0], [0], marker='o', color='w', markerfacecolor="#0072B2",
               markeredgecolor='black', markersize=9, label='Exogenous'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor="#56B4E9",
               markeredgecolor='black', markersize=9, label='Flexibility feature'),
        Line2D([0], [0], marker='o', color='w', markerfacecolor="#D55E00",
               markeredgecolor='black', markersize=9, label='Outcome'),
        Line2D([0], [0], color="#555555", lw=2, alpha=0.75,
               label='Undirected edge (thickness ∝ bootstrap freq.)'),
    ]
    ax.legend(handles=handles, loc='lower left', fontsize=7.5, frameon=True)
    ax.set_title(title, fontsize=11, loc='left')
    ax.set_axis_off()


def draw_stability_bars(ax, edges: pd.DataFrame, top_k: int = 12,
                        title: str = "(b) Bootstrap stability for the same undirected edges"):
    edges = edges.sort_values("frequency", ascending=False).head(top_k).copy()
    # Use '—' (em dash) as the explicit undirected separator
    labels = [f"{short(r['src'])}  —  {short(r['dst'])}" for _, r in edges.iterrows()]
    freqs = edges["frequency"].astype(float).tolist()

    y = np.arange(len(labels))[::-1]
    bars = ax.barh(y, freqs, color="#E69F00", edgecolor="black", lw=0.4)
    for yi, v in zip(y, freqs):
        ax.text(v + 0.01, yi, f"{v:.2f}", va='center', fontsize=7.5)

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("Bootstrap stability frequency (100 iterations)")
    ax.set_xlim(0, 1.08)
    ax.set_title(title, fontsize=11, loc='left')
    ax.grid(axis='x', alpha=0.3, linestyle='--')
    ax.set_axisbelow(True)


def render_consistent(edges: pd.DataFrame) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6.2),
                             gridspec_kw={"width_ratios": [1.05, 1.0]})
    draw_undirected_graph(axes[0], edges)
    draw_stability_bars(axes[1], edges)
    fig.suptitle(
        "Fig 6 (re-rendered, workstream #E). Both panels share the SAME underlying PC-bootstrap output "
        "(undirected edges). Arrowheads removed; '—' used to denote undirected adjacency in panel (b).",
        fontsize=9.5, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_BASE / "fig6_consistent.png", dpi=200, bbox_inches="tight")
    fig.savefig(OUT_BASE / "fig6_consistent.pdf", bbox_inches="tight")
    fig.savefig(OUT_BASE / "fig6_consistent.svg", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote fig6_consistent.{{png,pdf,svg}}")


def render_before_after(edges: pd.DataFrame) -> None:
    fig = plt.figure(figsize=(15, 8.5))
    gs = fig.add_gridspec(2, 2, hspace=0.35, wspace=0.15)

    # Row 1: published (before)
    ax_bef_a = fig.add_subplot(gs[0, 0])
    ax_bef_b = fig.add_subplot(gs[0, 1])
    if PUBLISHED_6A.exists():
        ax_bef_a.imshow(mpimg.imread(str(PUBLISHED_6A)))
    ax_bef_a.set_title("(BEFORE — published Fig 6a)\nArrowheads imply orientation",
                       fontsize=9.5, loc='left')
    ax_bef_a.set_axis_off()
    if PUBLISHED_6B.exists():
        ax_bef_b.imshow(mpimg.imread(str(PUBLISHED_6B)))
    ax_bef_b.set_title("(BEFORE — published Fig 6b)\n'→' labels imply direction",
                       fontsize=9.5, loc='left')
    ax_bef_b.set_axis_off()

    # Row 2: re-rendered (after)
    ax_aft_a = fig.add_subplot(gs[1, 0])
    ax_aft_b = fig.add_subplot(gs[1, 1])
    draw_undirected_graph(ax_aft_a, edges,
                          title="(AFTER) Undirected adjacency (PC output)")
    draw_stability_bars(ax_aft_b, edges,
                        title="(AFTER) Stability with '—' separator")

    fig.suptitle(
        "Fig 6 before / after — arrow-convention audit (workstream #E).  "
        "Underlying PC bootstrap emitted 100 % undirected edges; "
        "the re-rendered version (bottom row) removes arrowheads and replaces "
        "'→' with '—' to match algorithmic output.",
        fontsize=10.5, y=0.99)
    fig.savefig(OUT_BASE / "fig6_before_after.png", dpi=180, bbox_inches="tight")
    fig.savefig(OUT_BASE / "fig6_before_after.pdf", bbox_inches="tight")
    plt.close(fig)
    print(f"  wrote fig6_before_after.{{png,pdf}}")


def write_summary(edges: pd.DataFrame) -> None:
    # Verify undirected status of the cached PC output
    dag = pd.read_csv(DAG_EDGES_CSV)
    kinds = dict(dag["kind"].value_counts())
    n_total = int(len(dag))

    md = [
        "# Fig 6 arrow-consistency audit — summary",
        "",
        "## Underlying PC-bootstrap output (audit of cached results)",
        "",
        f"- `{DAG_EDGES_CSV.relative_to(u.ROOT)}` rows: {n_total}",
        f"- Edge kinds: {kinds}",
        ("- Conclusion: **all bootstrap-stable edges are undirected** (no orientation "
         "rule succeeded within the bootstrap-stable set produced by the PC algorithm "
         "with the fisherz independence test, alpha = 0.05, 100 bootstrap iterations)."),
        "",
        "## What Reviewer 2 flagged",
        "",
        ("> \"Why are arrows shown in Fig. 6b (L626) but not in Fig. 6a?\""),
        "",
        ("In the published Fig 6, Panel (a) renders arrowheads on the graph edges "
         "and Panel (b) labels each edge using the '→' notation. Both visual "
         "conventions imply orientation that the underlying PC-bootstrap output "
         "did not actually infer (every edge in `Fig05_causal_dag__dag_edges.csv` "
         "has `kind = undirected`)."),
        "",
        "## Fix (this workstream)",
        "",
        ("Both panels are re-rendered with an explicit **undirected convention**: "
         "no arrowheads in Panel (a) and the em-dash '—' separator in Panel (b). "
         "Edge thickness in Panel (a) encodes the same bootstrap frequency that "
         "Panel (b) reports as a horizontal-bar value, making the two panels "
         "visually consistent and faithful to the algorithmic output. Files:"),
        "",
        "- `fig6_consistent.{png,pdf,svg}` — clean two-panel side-by-side, undirected throughout",
        "- `fig6_before_after.{png,pdf}` — published Fig 6a/6b above, re-rendered consistent version below",
        "",
        "## Recommended caption (replaces the published Fig 6 caption)",
        "",
        ("See `docs/revision/captions/fig6_legend_v2.md` for the paste-ready caption draft. "
         "Briefly: Panel (a) is the *undirected adjacency network* from PC bootstrap with "
         "edges sorted by bootstrap stability frequency; Panel (b) is the same edge set "
         "shown as a stability bar chart. The PC algorithm produced no oriented edges in "
         "the bootstrap-stable set, so both panels use the '—' undirected convention."),
        "",
        "## Caveats (transparency)",
        "",
        ("- The two panels (top-K = 12) match each other but show a *subset* of the full "
         "edge_stability table (`Fig05_causal_dag__edge_stability.csv`, 21+ edges). The "
         "top-K cutoff is for readability, not a discovery step."),
        ("- The re-rendered figure is intended for **SI**, not as a substitute for the "
         "published main-text Fig 6 image. The Fig 6 caption in the manuscript will be "
         "updated according to `fig6_legend_v2.md` to reconcile the rendered figures "
         "with the algorithmic output."),
    ]
    (OUT_BASE / "fig6_audit_summary.md").write_text(
        "\n".join(md) + "\n", encoding="utf-8")
    print(f"  wrote fig6_audit_summary.md")


def main() -> None:
    print("[09 fig6 audit] starting...")
    edges = load_data()
    print(f"  edge_stability rows: {len(edges)}")
    render_consistent(edges)
    render_before_after(edges)
    write_summary(edges)
    print("[09 fig6 audit] done.")


if __name__ == "__main__":
    main()
