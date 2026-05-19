# Claude Code Handoff — iScience revision (2026-05-19)

This file is the **single source of truth** for what's done, what's
pending, and exactly what to type into Claude Code to continue.

## How to start

```powershell
# In Windows PowerShell:
cd C:\Cursor\Stenotrophomonas-causal AI
claude --permission-mode acceptEdits
```

When the session opens, **first read these three files in order** so the
agent has full context:

1. `CLAUDE.md`                              ← project constitution, hard constraints
2. `docs/revision/HANDOFF.md`               ← this file (state + pending work)
3. `revision_runs/iscience_rev1/REPORT.md`  ← already-completed workstream summaries

You can paste the entire prompt below to bring Claude Code up to speed:

```
Read CLAUDE.md, docs/revision/HANDOFF.md, and revision_runs/iscience_rev1/REPORT.md
in that order. Then list the pending tasks from HANDOFF.md ordered by priority
and ask me which one to start with. Do not modify any files yet.
```

---

## Branch & git state

- Active branch: `revision/iscience-rev1` (from `main`, 7+ commits ahead)
- Hard constraint: do NOT touch the 468 pre-existing uncommitted changes
  under `acetate_xai/` — those were intentional and the previous Claude
  Code session left them as-is per user instruction.
- Workstream #7 changes (this handoff) are uncommitted because the
  sandbox could not delete a stale `.git/index.lock`. Commit them first
  on Windows (see §A below).

---

## ✅ Completed workstreams

| # | Workstream | Headline | Files |
|---|---|---|---|
| 1 | Feature panel ablation | curated 42 vs 300 vs random-30: macro-F1 0.957–0.991, R² 0.90–0.92 — panel size insensitive | `revision_runs/iscience_rev1/01_feature_panel_ablation/` |
| 2 | Baseline benchmarking | inputs/GEM-summary/curated × LogReg/RF/XGB. Curated+RF F1=1.0, Curated+XGB F1=0.957 | `revision_runs/iscience_rev1/02_benchmarking/` |
| 3 | Existing-data perf summary | 5-fold CV macro-F1=0.991, R²=0.906; C10 = top mismatch | `revision_runs/iscience_rev1/05_existing_data/` |
| 4 | iML1515 transfer | macro-F1=0.972, R²=0.978 on 244 LHS conditions | `revision_runs/iscience_rev1/03_transfer/` |
| 5 | Runtime profiling | total ≈ 3.2 s, <420 MB peak RSS | `revision_runs/iscience_rev1/04_runtime/` |
| 7 | **Point-flux vs flexibility (novelty defense)** | iSO1: W 0.957 / M 0.941 / PF 0.939 / PFA 0.972 / OBJ 0.912. iML1515: W 0.973 / OBJ 0.289. Honest: W beats same-FVA midpoint by +0.016 F1; novelty is **conceptual**, not raw F1 dominance | `revision_runs/iscience_rev1/07_pointflux/` |

All rebuttal-ready Results / Methods / Editor-response paragraphs are
already drafted in `docs/revision/rebuttal_insertions.md`.

---

## 🟡 Pending — in priority order

### A. **Commit workstream #7** (1 min — do this first)

```powershell
cd C:\Cursor\Stenotrophomonas-causal AI
git add code/revision/07_pointflux_baseline.py code/revision/_pfba_runner.py `
        revision_runs/iscience_rev1/07_pointflux/ `
        revision_runs/iscience_rev1/figures/pointflux_vs_width.* `
        docs/revision/rebuttal_insertions.md `
        revision_runs/iscience_rev1/REPORT.md `
        revision_runs/iscience_rev1/metrics_summary.csv `
        docs/revision/HANDOFF.md
git commit -m "revision: point-flux vs flexibility-interval baseline (novelty defense)"
```

### B. **Fill in iML1515 pFBA baseline** (~3 min — only needs the user's already-cached iML1515)

The sandbox could not download iML1515 from BiGG. The Windows side
already has it cached. Run:

```powershell
cd C:\Cursor\Stenotrophomonas-causal AI
python code/revision/_pfba_runner.py iml 0 244
python code/revision/07_pointflux_baseline.py
```

Expected: `pfba_fluxes_iml1515.parquet` appears in `07_pointflux/`, and
the iML1515 columns of `pointflux_vs_width.png` are now filled.
Re-commit with: `git commit -am "revision: iML1515 pFBA baseline (completes workstream #7)"`.

### C. **Conceptual Table 1 / Fig 1A** (≈ 30 min in Claude Code)

Goal: a one-glance visual that contrasts **point-state prediction** with
**feasible-space diagnosis**. Place at the start of the manuscript so
Reviewer 1's "I can't see the framework logic from the abstract" is
solved on page 1.

**Paste this prompt into Claude Code:**

```
Create a conceptual Table 1 OR Figure 1A (your choice; whichever reads
faster) that contrasts the two diagnostic paradigms head-to-head.

Required columns / rows:
- Diagnostic question (point vs feasible-space)
- Mathematical object (one flux vector vs interval [vmin, vmax])
- What is interpreted (magnitude vs width)
- Mechanistic reading (optimum location vs degrees of freedom / rigidification)
- Downstream artefact (single flux map vs rigidification/causal map)
- Examples (pFBA, dFBA vs ours)

Constraints:
- Build under code/revision/ and write the figure to
  revision_runs/iscience_rev1/figures/conceptual_table1.{png,pdf,svg}
- If you make it a figure, also save an editable SVG so the user can
  paste into the Word manuscript.
- Do NOT modify the manuscript .docx — only create the figure/table.
- Keep ≤ 6 rows × 3 columns; use a clean black-and-white style suitable
  for iScience.
- Save a short markdown caption draft at
  docs/revision/captions/conceptual_table1.md
```

### D. **Fig 5 non-linearity reanalysis** (≈ 1 hr)

Reviewer 2: "Fig 5b/c look linear even though you call them non-linear."

**Paste this prompt:**

```
Reviewer 2 questioned the "non-linearity" claim around Fig 5. We need a
defense WITHOUT touching the manuscript figures themselves.

Tasks:
1. Recompute SHAP dependence (Booster.predict pred_contribs=True) for
   the same severity-regression model used in the paper Fig 5 panels.
   Reuse acetate_xai/scripts/train_xgb_shap_severity.py outputs if
   convenient, otherwise reproduce.
2. Overlay LOESS / lowess fit (statsmodels) on each Fig-5-equivalent
   dependence plot. Also fit a 2-knot piecewise-linear regression and
   report the F-test against a single linear fit.
3. For the top-2 severity SHAP features, compute SHAP interaction
   values (XGBoost native) and visualize. A clear interaction or a
   piecewise break is what justifies "non-linear effects".
4. Write three deliverables:
   - revision_runs/iscience_rev1/08_nonlinearity/fig5_dependence_with_fits.{png,pdf}
   - revision_runs/iscience_rev1/08_nonlinearity/nonlinearity_summary.md
     (with F-test p-values and a one-sentence honest verdict)
   - A paste-ready paragraph appended to docs/revision/rebuttal_insertions.md
     under a new "## 8. Non-linearity defense (Reviewer 2)" section.

Conservative framing if the test fails: do NOT claim non-linearity in
the manuscript wording — instead reframe as "threshold-like or
interaction-modulated effects" and update Fig 5 caption language in
docs/revision/captions/fig5_caption_v2.md (do not edit the docx).
```

### E. **Fig 6 arrows audit** (≈ 20 min)

Reviewer 2: "Fig 6a has no arrows but Fig 6b does — clarify."

**Paste this prompt:**

```
Reviewer 2 flagged Fig 6a/6b inconsistency (arrows vs no arrows). Our
strategic answer: Fig 6a should be relabeled / explained as "undirected
adjacency" while Fig 6b is "directed edge stability". Both panels come
from the same PC-bootstrap procedure but represent different stages of
the causal discovery output.

Tasks:
1. Read acetate_xai/scripts/run_causal_discovery.py (or equivalent) and
   confirm the actual output: undirected skeleton (a) vs directed
   acyclic edges (b).
2. Generate an updated SVG legend / caption draft in
   docs/revision/captions/fig6_legend_v2.md that makes the distinction
   explicit.
3. Re-render Fig 6 panels with consistent arrow conventions to
   revision_runs/iscience_rev1/09_fig6_audit/fig6_consistent.{png,pdf,svg}.
   If you must reuse paper Fig 6 sources, save a side-by-side
   before/after to make the change traceable.
4. Append a "## 9. Fig 6 arrow consistency" section to
   docs/revision/rebuttal_insertions.md with a 2-sentence rebuttal.

Do NOT edit the manuscript figure files directly.
```

### F. **Point-by-point response letter** (≈ 45 min)

**Paste this prompt:**

```
Generate a complete point-by-point response letter at
docs/revision/response_letter.docx (use the docx skill).

Source material:
- Submission/Revision_1.pdf — the editor and reviewer comments
- docs/revision/rebuttal_insertions.md — all completed workstream answers
- revision_runs/iscience_rev1/REPORT.md — quantitative results

Structure required by iScience:
1. One-page cover paragraph (reposition novelty as feasible-space
   diagnosis vs point-state prediction; mention every new workstream
   in one sentence each)
2. Per-comment table with three columns:
   a) reviewer/editor quote (copy verbatim from Revision_1.pdf)
   b) our response (1-2 paragraphs, conservative wording, cite line
      numbers of the revised manuscript using placeholders like
      [Line XXX])
   c) corresponding revision file / figure / SI section
3. Final "summary of changes" table:
   - New figures added (list with placement: main vs SI)
   - New SI sections
   - References added
   - Caveats now stated in main text

Constraints:
- Use the docx skill at
  C:\Users\user\AppData\Roaming\Claude\local-agent-mode-sessions\skills-plugin\...\skills\docx\SKILL.md
  Read its instructions BEFORE building the doc.
- Conservative wording throughout: "proof-of-concept", "hypothesis-
  prioritization", "in silico demonstration" — do NOT claim
  "validated tool" or "deployable system".
- Distinguish between mandatory and recommended responses; mandatory
  ones must each have a concrete pointer to a revision artefact.
```

### G. **Novelty paragraphs for manuscript** (deferred until F is done)

These are the paste-ready blocks for the Word doc itself. They're
already drafted in `docs/revision/rebuttal_insertions.md` §7 — when the
user is ready to do the actual Word edits, just copy from there:

- Abstract last sentence (feasible-space framing)
- Introduction new paragraph (the one Reviewer 1 asked for)
- Discussion opener (3 sentences)
- Limitations paragraph (claim refinement per editor)

---

## 📐 Environment notes

- Python: `C:/Users/user/Miniforge3/python.exe` (3.12.2)
- Solver: optlang.glpk_interface (default)
- Critical packages: cobra 0.31.1, xgboost 3.1.3, scikit-learn 1.8.0,
  shap 0.49.1 (works after `pip install "numpy<2.4"`),
  causal-learn 0.1.4.5
- iML1515 is **already cached** on the Windows machine via
  `cobra.io.load_model('iML1515')` — no re-download needed.
- The sandbox-side runs in this revision used Python 3.10.12 +
  xgboost 3.2.0 in the Linux sandbox; results are consistent because
  the dataset replay matches the original `objective_value` exactly
  (verified Δ=0.0 across 8 random rows).

---

## 🛑 Hard rules for Claude Code (re-confirm before any edit)

From `CLAUDE.md`:

- Do NOT touch `_unpacked_main/`, `_unpacked_si/`, the published
  release archives, `BaseModel.xml`, `Final_genome.fasta`, or
  existing committed scripts under `acetate_xai/scripts/`.
- Add new revision scripts under `code/revision/` only.
- All outputs under `revision_runs/iscience_rev1/` only.
- Determinism: fixed `--seed 42` (random-control panels: seed grid 1..10).
- After each workstream, append a 3-bullet summary to `REPORT.md` and a
  paste-ready section to `docs/revision/rebuttal_insertions.md`.
- Do NOT modify the Word/PDF manuscript files directly. All manuscript
  edits go through the user, who will paste from the rebuttal/handoff
  drafts.

---

## ⏭ Suggested order for the next Claude Code session

1. A (commit #7) — 1 min, no thinking
2. B (iML1515 pFBA) — 3 min, runs and re-commits
3. F (response letter) — highest leverage; it ALSO surfaces what extra
   evidence each reviewer point still needs
4. C (conceptual Table 1) — fast and very high-impact (page-1 visibility)
5. D + E (Fig 5/6 fixes) — answer Reviewer 2's last two specific points

After F, the user has a complete submission package: rebuttal letter +
all SI figures + numbered revision artefacts. Manuscript-side text
edits (G) are then a paste job, not an analysis job.
