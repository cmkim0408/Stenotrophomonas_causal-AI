# Fig 6 arrow-consistency audit — summary

## Underlying PC-bootstrap output (audit of cached results)

- `results\figures_final\data\Fig05_causal_dag__dag_edges.csv` rows: 12
- Edge kinds: {'undirected': np.int64(12)}
- Conclusion: **all bootstrap-stable edges are undirected** (no orientation rule succeeded within the bootstrap-stable set produced by the PC algorithm with the fisherz independence test, alpha = 0.05, 100 bootstrap iterations).

## What Reviewer 2 flagged

> "Why are arrows shown in Fig. 6b (L626) but not in Fig. 6a?"

In the published Fig 6, Panel (a) renders arrowheads on the graph edges and Panel (b) labels each edge using the '→' notation. Both visual conventions imply orientation that the underlying PC-bootstrap output did not actually infer (every edge in `Fig05_causal_dag__dag_edges.csv` has `kind = undirected`).

## Fix (this workstream)

Both panels are re-rendered with an explicit **undirected convention**: no arrowheads in Panel (a) and the em-dash '—' separator in Panel (b). Edge thickness in Panel (a) encodes the same bootstrap frequency that Panel (b) reports as a horizontal-bar value, making the two panels visually consistent and faithful to the algorithmic output. Files:

- `fig6_consistent.{png,pdf,svg}` — clean two-panel side-by-side, undirected throughout
- `fig6_before_after.{png,pdf}` — published Fig 6a/6b above, re-rendered consistent version below

## Recommended caption (replaces the published Fig 6 caption)

See `docs/revision/captions/fig6_legend_v2.md` for the paste-ready caption draft. Briefly: Panel (a) is the *undirected adjacency network* from PC bootstrap with edges sorted by bootstrap stability frequency; Panel (b) is the same edge set shown as a stability bar chart. The PC algorithm produced no oriented edges in the bootstrap-stable set, so both panels use the '—' undirected convention.

## Caveats (transparency)

- The two panels (top-K = 12) match each other but show a *subset* of the full edge_stability table (`Fig05_causal_dag__edge_stability.csv`, 21+ edges). The top-K cutoff is for readability, not a discovery step.
- The re-rendered figure is intended for **SI**, not as a substitute for the published main-text Fig 6 image. The Fig 6 caption in the manuscript will be updated according to `fig6_legend_v2.md` to reconcile the rendered figures with the algorithmic output.
