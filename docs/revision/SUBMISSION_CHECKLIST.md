# Submission readiness — checklist

**Manuscript:** ISCIENCE-D-26-04043
**Branch:** `revision/iscience-rev1`
**Date:** 2026-05-19

This is the **single user-facing checklist** for the revision submission.
Work through sections A → E in order. Every paste source is in
`docs/revision/SUBMISSION_BUNDLE/` (auto-generated) or in
`docs/revision/rebuttal_insertions.md` (the master draft file).

---

## At a glance

- **Quantitative analyses:** all 6 revision workstreams (#1, #2, #3, #4, #5, #7) plus the extended FVA campaign and the #D, #E follow-ups are complete. Deferred-count = **0**; every editor mandatory revision and every reviewer comment is answered in this iteration.
- **Auto-generated deliverables:** response letter, highlights, graphical-abstract brief, Conceptual Table 1 spec, AI disclosure, STAR Methods guide, cover-letter paragraph, caption replacements — all in `SUBMISSION_BUNDLE/`.
- **Housekeeping verify:** 6/7 items already present in the main draft; 1 item(s) need a small addition (see §E).
- **Estimated user time after this checklist is generated:** ~1.5 hours (paste edits ~30 min + cover letter ~5 min + graphical abstract ~30 min + EM upload ~20 min).

---

## A. Manuscript text edits (paste from `rebuttal_insertions.md`)

Open the Word manuscript with **Track Changes ON**. Paste each block at
the placeholder line indicated. The `[Line XXX]` markers should be
replaced with the actual line numbers in your Word doc as you paste
(the rebuttal letter references `[Line XXX]` placeholders the editor
will accept; you don't need them to be real until you finalise the
highlighted version in §B).

- [ ] **A1.** Abstract — append 1–2 sentences emphasising the
  feasible-space framing. Source: `rebuttal_insertions.md` §7
  Cover-letter insert (first paragraph).
- [ ] **A2.** Introduction — insert the **new paragraph** after current
  paragraph 2. Source: `rebuttal_insertions.md` §7 Introduction insert.
- [ ] **A3.** Methods (Feature engineering via targeted FVA) — add the
  extended-FVA campaign paragraph. Source: `rebuttal_insertions.md` §1
  Methods paragraph + the *Width universe — extended FVA campaign*
  bullet from `REPORT.md` Known limitations.
- [ ] **A4.** Methods (Explainable AI and causal discovery) — add the
  **feature-panel ablation** paragraph (`rebuttal_insertions.md` §1
  Methods) and the **baseline benchmarking** paragraph (§2 Methods).
- [ ] **A5.** Methods (QUANTIFICATION AND STATISTICAL ANALYSIS) — add
  the runtime / environment 2–3 sentences. Source: §5 Methods.
- [ ] **A6.** Methods (Explainable AI and causal discovery) — add the
  **point-flux vs flexibility-interval** paragraph. Source: §7 Methods.
- [ ] **A7.** Methods — add the **Fig 5 non-linearity reanalysis**
  paragraph. Source: §8 Methods.
- [ ] **A8.** Results (near Fig 7) — add the existing-data 5-fold CV
  sentences (macro-F1 = 0.991; R² = 0.906; C10 top mismatch).
  Source: §3 Results sentences.
- [ ] **A9.** Discussion — add the transfer / generalizability
  paragraph. Source: §4 Results paragraph.
- [ ] **A10.** Discussion — add the **3-sentence opener** that
  re-positions the novelty as feasible-space diagnosis. Source: §7
  Results paragraph + Cover-letter insert.
- [ ] **A11.** Limitations — add the proof-of-concept / hypothesis-
  prioritization framing. Source: §7 Cover-letter insert (paragraph 3)
  + Strategic framing.
- [ ] **A12.** **Fig 5 caption** — replace with Variant A from
  `SUBMISSION_BUNDLE/caption_replacements/fig5_caption_v2.md`
  (STRONG_NONLINEAR verdict; explicit breakpoints + ΔR²).
- [ ] **A13.** **Fig 6 caption** — replace with Variant B from
  `SUBMISSION_BUNDLE/caption_replacements/fig6_legend_v2.md`
  (undirected convention; conditional-dependency map).

---

## B. Build the highlighted + clean Word files

- [ ] **B1.** After completing §A with Track Changes ON, save the file as
  `Main_Stenotrophomonas Causal AI_Main draft_REVISED-r2-highlighted.docx`.
- [ ] **B2.** Duplicate the highlighted file, accept all tracked
  changes, save as `Main_..._REVISED-r2-clean.docx`.

iScience requires **both** versions on submission.

---

## C. Update the cover letter

- [ ] **C1.** Open existing
  `Submission/Causal AI-Cover letter-final-iScience.docx`.
- [ ] **C2.** Insert the `cover_letter_paragraph.md` text
  (`SUBMISSION_BUNDLE/cover_letter_paragraph.md`) as the **second**
  paragraph of the cover (after the brief restatement of the
  manuscript scope, before the closing reviewer acknowledgement).
- [ ] **C3.** Save as `Causal AI-Cover letter-r2.docx`.

---

## D. New deliverables (auto-spec'd; user PowerPoint)

- [ ] **D1. Highlights (3–5 bullets).** Use the compressed `<85`-char
  variants from `SUBMISSION_BUNDLE/highlights_bullets.md`. Paste into
  the iScience EM submission form *and* keep a copy near the Abstract.
- [ ] **D2. Graphical abstract.** Follow the three-panel design brief
  in `SUBMISSION_BUNDLE/graphical_abstract_brief.md` to draw a PNG
  in PowerPoint (~30 min). Export at 600 dpi.
- [ ] **D3. (Optional) Conceptual Table 1 / Fig 1A.** Spec at
  `SUBMISSION_BUNDLE/conceptual_table1_spec.md`. Reviewer 1's intro-
  clarity comment is already addressed narratively in §A2; Table 1
  is a visual bonus.
- [ ] **D4. AI disclosure (optional refinement).** Existing manuscript
  paragraph already discloses ChatGPT use. If you want to also disclose
  Claude Code usage during the revision, swap in the refined paragraph
  from `SUBMISSION_BUNDLE/ai_disclosure_paragraph.md`.

---

## E. Editor housekeeping verify

Auto-verified status (6/7 present). Full
evidence table in `SUBMISSION_BUNDLE/housekeeping_verify.md`.

- [x] **H_DATA** — Data and code availability statement — present in main draft (see housekeeping_verify.md)
- [x] **H_CREDIT** — CRediT author statement (Author contributions) — present in main draft (see housekeeping_verify.md)
- [x] **H_COI** — Conflict of interest / Declaration of interests — present in main draft (see housekeeping_verify.md)
- [x] **H_FUNDING** — Funding statement — present in main draft (see housekeeping_verify.md)
- [ ] **H_ETHICS** — Ethics approval (if applicable) — No human or animal subjects involved (bacterial culture study). iScience editor's instruction: 'Ethics approval (if applicable)' — recommendation: add a single explicit sentence under Resource availability noting that no ethics approval is required as the study uses only a bacterial isolate.
- [x] **H_AI** — Generative AI / AI-assisted tools disclosure — present in main draft (see housekeeping_verify.md)
- [x] **H_STAR** — STAR Methods format — present in main draft (see housekeeping_verify.md)

---

## F. Final upload (editorialmanager.com/iscience)

Files to upload (zip `SUBMISSION_BUNDLE/` is not directly accepted
by EM; you'll upload individual files):

- [ ] **F1.** Cover letter — `Causal AI-Cover letter-r2.docx` (§C3).
- [ ] **F2.** Manuscript (highlighted) — `Main_..._REVISED-r2-highlighted.docx` (§B1).
- [ ] **F3.** Manuscript (clean) — `Main_..._REVISED-r2-clean.docx` (§B2).
- [ ] **F4.** Response letter — `SUBMISSION_BUNDLE/response_letter.docx`.
- [ ] **F5.** Supplementary Information — your existing SI .docx +
  the new SI figures from `revision_runs/iscience_rev1/figures/` and
  `revision_runs/iscience_rev1/0[1-9]_*/`. Include:
    - `feature_panel_ablation.{png,pdf}`
    - `benchmark_comparison.{png,pdf}`
    - `confusion_matrix.{png,pdf}`
    - `regression_residuals.{png,pdf}`
    - `external_transfer.{png,pdf}`
    - `pointflux_vs_width.{png,pdf}`
    - `08_nonlinearity/fig5_dependence_with_fits.{png,pdf}`
    - `08_nonlinearity/fig5_shap_interactions.{png,pdf}`
    - `09_fig6_audit/fig6_consistent.{png,pdf,svg}`
    - `09_fig6_audit/fig6_before_after.{png,pdf}`
  Plus the SI tables (`*_metrics.csv`, `runtime_summary.csv`, etc.).
- [ ] **F6.** Graphical abstract — `graphical_abstract.png` (§D2).
- [ ] **F7.** Highlights — paste into the EM submission form (§D1).
- [ ] **F8.** Verify EM author list matches submission (no changes).

---

## Ready-to-upload status

**Auto-generated deliverables:** ✅ complete (7 housekeeping items checked, 6 present, 1 need user attention).

**Remaining user-side work** (~1.5 h): manuscript paste (§A; ~30 min) → highlighted + clean Word files (§B; ~5 min) → cover letter update (§C; ~5 min) → graphical abstract image (§D2; ~30 min) → EM upload (§F; ~20 min).

After §A → §F complete, the submission is **READY TO UPLOAD** to `https://www.editorialmanager.com/iscience/`.
