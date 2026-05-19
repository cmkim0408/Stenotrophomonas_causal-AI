# Fig 6 — revised caption (workstream #E)

This file contains a paste-ready replacement caption for Fig 6 in the revised
main text. Two variants are provided: a **conservative SI-only version**
(recommended) and an **in-place main-text replacement** for the case where
the user prefers to update the published Fig 6 caption directly.

The change is *caption-level only*; the underlying PC-bootstrap data are
unchanged. Both variants explicitly state that the bootstrap-stable edge set
is undirected, reconciling the visual representation with the algorithmic
output (`Fig05_causal_dag__dag_edges.csv`: 100 % of rows have
`kind = undirected`).

---

## Variant A — SI caption (recommended)

> **Figure SX. Re-rendered causal-discovery panels (undirected convention; SI to Fig 6).**
> (a) **Undirected adjacency network** inferred by the PC algorithm with the
> Fisher-z independence test (α = 0.05) under 100 bootstrap resamples. Each
> edge connects a pair of variables whose conditional dependence survived the
> bootstrap procedure with frequency ≥ 0.10; edge thickness encodes the
> bootstrap stability frequency. **Arrowheads have been omitted**: the
> bootstrap-stable edge set produced by PC contained no edges that the
> orientation rules could direct (the equivalence class is fully undirected
> for this dataset).
> (b) **Bootstrap stability** for the same edge set, sorted by frequency.
> Edge labels use the em-dash ("—") separator to denote the undirected
> adjacency; the value at the right of each bar is the fraction of 100
> bootstrap resamples in which that edge survived.
> Variable colour coding in (a): blue = exogenous (ATPM-fixed, O₂ uptake
> bound); cyan = SHAP flexibility feature; orange = outcome (regime,
> growth severity). The same top-12 edges are shown in both panels; the full
> stability table is provided as Table S6.

---

## Variant B — In-place main-text Fig 6 caption replacement

> **Figure 6. Inferred rigidification structure from PC-bootstrap causal discovery.**
> (a) **Adjacency graph** of the top-K conditional-dependency edges that
> survived 100 bootstrap resamples of the PC algorithm (Fisher-z
> independence test, α = 0.05; stability cutoff = 0.10). Edges are rendered
> as **undirected** because the bootstrap-stable equivalence class produced
> by PC for this dataset contained no oriented edges; edge thickness encodes
> bootstrap stability frequency. Variable colour coding: blue = exogenous
> (ATPM-fixed, O₂ uptake bound); cyan = SHAP flexibility feature; orange =
> outcome (regime, growth severity).
> (b) **Bootstrap stability** for the same edge set, expressed as the fraction
> of 100 resamples in which the edge survived. Edge labels use the em-dash
> ("—") to indicate undirected adjacency. The same top-K edges shown in (a)
> are reproduced in (b) so that the two panels are visually consistent.
> Together the panels show that the inferred rigidification structure is a
> conditional-dependency map rather than a fully oriented causal DAG, in
> agreement with the manuscript's interpretation of the rigidification map
> as a **hypothesis-prioritization structure**, not causal proof.

---

## Note for the response letter

The two-paragraph rebuttal in §9 of `docs/revision/rebuttal_insertions.md`
explicitly describes the audit (cached PC output kinds, re-render with
undirected convention, no algorithmic re-run) and cites
`revision_runs/iscience_rev1/09_fig6_audit/fig6_consistent.{png,pdf,svg}`
and `fig6_before_after.{png,pdf}`.

The same response also closes editor mandatory revision **M7** and reviewer
2's comment **R2.2**.
