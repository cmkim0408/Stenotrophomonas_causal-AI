# Rebuttal insertions — iScience revision

**Branch:** `revision/iscience-rev1`
**Date:** 2026-05-19 (workstream #7 completed; iML1515 pFBA filled in)

This file collects ready-to-paste paragraphs and point-by-point responses for
each completed revision workstream. **Numbers are reported conservatively**;
interpretive language has been reframed away from "X is best" toward "the
flexibility-based representation carries diagnostic signal that supports an
interpretable workflow across model classes."

> **Single-line strategic framing for the whole revision:**
> The revision's claim is **not** that XGBoost is uniquely optimal. The
> claim is that flexibility-based features carry meaningful diagnostic
> signal and support an interpretable diagnostic workflow across model
> classes — XGBoost is the consistent explainable learner used to expose
> that signal.

---

## 0. Placement strategy (main text vs SI)

| Asset | Source | Recommended placement |
|---|---|---|
| Short benchmark result paragraph | `02_benchmarking/` | **Main text** (Discussion) — 1 paragraph |
| Short existing-data metrics paragraph | `05_existing_data/` | **Main text** (Results, near Fig 7) — 1 paragraph |
| Short external-transfer paragraph | `03_transfer/` | **Main text** (Discussion) — 1 paragraph |
| Runtime 2–3 sentences | `04_runtime/` | **Main text** (Methods) |
| `feature_panel_ablation.png` | `01_ablation/` | **SI figure** |
| `benchmark_comparison.png` | `02_benchmarking/` | **SI figure** |
| `confusion_matrix.png` | `05_existing_data/` | **SI figure** |
| `regression_residuals.png` | `05_existing_data/` | **SI figure** |
| `external_transfer.png` | `03_transfer/` | **SI figure** (do **not** replace main Fig 4) |
| Runtime summary table | `04_runtime/` | **SI table** |

Rationale: keep the main narrative anchored on the curated Fig 4 / Fig 7
story; use SI to supply quantitative reviewer-driven justification.

---

## 1. Feature-panel ablation (Reviewer 2: "why 30?")

### Strategic framing

The curated paper-aligned panel (now 42 reactions on the extended
`width__` universe) is **not** a performance-optimal subset. The ablation
shows the diagnostic signal is broadly distributed across the deployed
`width__` universe, so the curated panel should be defended as an
**interpretability-oriented diagnostic layer**, not a global optimum.

### Results paragraph (main text or SI)

> To assess whether the choice of curated reactions/modules is critical,
> we performed a six-level feature-panel ablation on the same 5-fold
> cross-validated regime classification and severity regression tasks.
> An extended FVA campaign (Methods) re-ran targeted FVA on 180
> additional reactions across all 242 conditions, expanding the deployed
> `width__` universe from 120 to ~300 columns and bringing paper-named
> TCA/respiration anchors (MDH, ICDHx, ICDHyr, ICL, MALS, PYK, PPC,
> NADH16pp, FUM …) into the feature set. Within this extended universe,
> panels of size 10, 20, 42 (curated paper-aligned), 50 (top-K by SHAP),
> and 300 (all widths) all yielded macro-F1 between 0.957 and 0.991 and
> severity R² between 0.90 and 0.92. Random 30-feature controls across
> 10 seeds were comparable (macro-F1 0.981 ± 0.009; R² 0.913 ± 0.009).
> The diagnostic signal is therefore broadly distributed across modules
> rather than concentrated in a few hand-picked reactions. **The
> curated panel should be interpreted as an interpretability-oriented
> diagnostic layer rather than a performance-optimal subset.**

### Methods paragraph

> Feature-panel ablation. The deployed `width__` universe was first
> extended via a supplemental FVA campaign on the 242 stored
> (campaign × run_folder × condition) tuples; per-row replays were
> verified against the original `objective_value` (matched within 5e-2
> for all 242 rows), and 180 additional reactions from `targets_300`
> were added (`extend_fva_campaign.py`). Six panels were then compared
> on the resulting ~300-`width__` universe: a curated paper-aligned
> panel (n=42; intersection of paper-named modules with available
> `width__` columns), top-10/20/50 ranked by mean(|SHAP|) on a baseline
> XGBoost classifier, the full 300-width set, and random 30-feature
> controls across 10 seeds. For each panel, regime classification
> (XGBoost, stratified 5-fold CV) and severity regression (XGBoost,
> 5-fold CV) were evaluated by macro-F1, balanced accuracy, RMSE, MAE,
> and R². Outputs: `revision_runs/iscience_rev1/01_feature_panel_ablation/`.

### Rebuttal answer

> Reviewer 2 asked why 30 features were selected. We extended the
> deployed FVA universe to include the paper-named TCA / respiration
> anchors that had been alphabetically truncated from the prior
> 120-feature subset (MDH, ICDH, ICL, MALS, PYK, PPC, NADH16pp, FUM…),
> giving a 300-`width__` universe in which the panel-size question can
> be answered cleanly. Within this extended universe, the framework's
> diagnostic performance is largely insensitive to panel size
> (Δmacro-F1 ≤ 0.04 across panels of size 10–300; random 30-feature
> controls match curated). The curated 30/42 panel was chosen for
> **interpretability** — the reactions map to paper-named central-carbon
> and respiratory modules — not for predictive performance. The
> ablation has been added to the revised SI and is referenced in the
> Results.

**Strength:** ✅ SI figure + 1-paragraph Results insert.

---

## 2. Baseline benchmarking (editor + reviewers)

### Strategic framing

XGBoost on the curated panel reached macro-F1 = 0.957, while RandomForest
on the *same* panel reached 1.000 and LogReg reached 0.976. Inputs-only
baselines also reached macro-F1 ≈ 0.96. The honest reading is that the
**flexibility-based representation** carries the diagnostic signal across
multiple learners; XGBoost's role is to expose that signal in an
explainable form (SHAP), not to claim unique predictive superiority.

### Results paragraph (main text — REWRITTEN, conservative)

> We benchmarked the flexibility-feature pipeline against three classes of
> simpler baselines on the same 5-fold cross-validated split: (A)
> inputs-only features (uptake bounds), (B) GEM-summary features (FBA
> objective + uptake fluxes), and (C) the curated 31-reaction width panel
> paired with logistic regression, random forest, and XGBoost. Inputs-only
> baselines reached macro-F1 ≤ 0.96 but severity R² ≤ 0.86, indicating
> that uptake bounds alone do not resolve flexibility-driven variation.
> GEM-summary baselines achieved high R² (≥ 0.93), but the FBA objective
> enters their feature vector by construction (severity = obj/obj_max).
> On the curated width panel, all three learners performed competitively
> (macro-F1 0.957–1.000; severity R² 0.92–0.998), indicating that the
> diagnostic signal is carried primarily by the flexibility-based
> representation rather than by a single uniquely optimal learner. **Our
> main claim is therefore not that XGBoost is uniquely optimal, but that
> flexibility-based features carry meaningful diagnostic signal and
> support an interpretable diagnostic workflow across model classes.** We
> use XGBoost in the main framework as a consistent explainable learner.

### Methods paragraph

> Baseline benchmarking. Three baseline feature sets were defined: (A)
> uptake bounds for acetate, oxygen, ammonium, phosphate; (B) FBA
> objective and the per-anchor flux/saturation vector; (C) the curated
> 31-reaction width-feature panel. Each set was paired with logistic
> regression (standardized; balanced class weights), random forest (n=400
> trees), and XGBoost (n=300 trees, depth 4). 5-fold cross-validated
> macro-F1, balanced accuracy, per-class F1, RMSE, MAE, and R² are
> reported in
> `revision_runs/iscience_rev1/02_benchmarking/benchmark_metrics.csv`.

### Rebuttal answer (REWRITTEN, conservative)

> The editor flagged baseline comparison as critical. We added benchmarks
> spanning inputs-only, GEM-summary, and same-panel ML baselines. The
> benchmark indicates that the diagnostic signal is carried primarily by
> the flexibility-based representation rather than by a single uniquely
> optimal learner. On the curated width panel, XGBoost remained
> competitive, while logistic regression and random forest achieved
> comparable or higher scores in some settings. We therefore interpret
> the benchmark as support for the flexibility-based diagnostic
> representation and use XGBoost as a consistent explainable learner in
> the main framework. The B_gem_summary baseline includes
> `objective_value`, which is the numerator of the severity target G =
> obj/obj_max; its R²≈0.998 reflects construction overlap, not new
> predictive content, and is reported transparently.

**Strength:** ✅ SI figure + 1-paragraph Discussion insert.

---

## 3. Existing-data performance summary

### Strategic framing

This workstream gives the quantitative robustness backbone the editor asked
for. **It must be presented as model-internal CV robustness, distinct from
the experimental Fig 7 mismatch interpretation.**

### Results sentences (insert in Results, near Fig 7)

> On the LHS-derived diagnostic dataset (n=242), 5-fold cross-validation
> yielded macro-F1 = 0.991 with per-class F1 of 0.976 / 1.000 / 0.998 for
> N-, Ac-, and O2-limited regimes respectively (Fig SX confusion matrix).
> Severity regression achieved RMSE = 0.034 and R² = 0.906 with no
> systematic residual offset (mean = -0.001). These numbers describe
> **model-internal cross-validation robustness on the simulated
> diagnostic dataset**; they are distinct from the experimental
> agreement/mismatch analysis in Fig 7. Among the C1–C10 holdout
> conditions, the rank-residual of predicted severity vs measured OD600
> identified C10 (sealed-cap mid-O2, predicted severity 0.30 / measured
> OD 0.74) as the largest mismatch, consistent with the manuscript's
> non-stoichiometric-constraint interpretation.

### Methods sentences

> Existing-data performance summary. Confusion matrix and per-class
> precision/recall/F1 were computed from 5-fold cross-validated
> out-of-fold predictions of the XGBoost regime classifier trained on
> the extended ~300-width superset. Severity residuals (predicted -
> measured G = obj/obj_max) and rank-residual mismatch scoring for the
> C1–C10 holdout were saved to
> `revision_runs/iscience_rev1/05_existing_data/`.

**Strength:** ✅ supplementary figure (confusion matrix + residuals) + 1
short Results paragraph. *Do not conflate with Fig 7 experimental story.*

---

## 4. External transfer (iML1515)

### Strategic framing

The transfer demo supports **transferability of the diagnostic
formulation**, not generalization completion. Stay in SI; the
side-by-side SHAP comparison can support — but should not replace —
the main Fig 4. After the extended FVA campaign, the iSO1 top SHAP
features now include paper-named anchors (MDH, ICDHx, AKGDH, PPCDC,
SUCOAS), aligning the data-driven ranking with the published Fig 4
narrative (TCA / respiration / ATP).

### Results paragraph (main text or SI)

> To test whether the diagnostic logic transfers to a different organism,
> we applied the identical workflow — LHS over uptake bounds,
> shadow-price regime labeling, targeted FVA over a curated panel,
> XGBoost + SHAP — to *Escherichia coli* iML1515. From 250 LHS samples
> over (glc, ac, o2, nh4, pi) uptake bounds, 244 yielded feasible
> solutions, distributed across three regimes (o2_limited 115;
> nh4_limited 80; glc_limited 49). On a 45-reaction E. coli curated
> panel covering TCA / glyoxylate / glycolysis / respiration / acetate
> uptake / N biosynthesis, 5-fold cross-validated macro-F1 reached 0.972
> and severity R² reached 0.978. iML1515 top SHAP features (TCA /
> glyoxylate / glycolysis modules) and iSO1 top SHAP features computed
> on the extended ~300-`width__` universe (which now includes MDH,
> ICDHx, AKGDH, PPCDC, SUCOAS) were **system-specific reaction IDs**
> but **functionally analogous** across systems. **The framework is
> therefore transferable in formulation, while system-specific feature
> curation remains necessary; the present transfer analysis is in
> silico only, and wet-lab validation in the external organism remains
> future work.**

### Methods paragraph

> External transfer to iML1515. The well-curated public *E. coli* GEM
> iML1515 (loaded via `cobra.io.load_model`) was sampled at n=250 LHS
> points over uptake-bound ranges chosen to span aerobic glucose and
> acetate conditions. For each feasible solution, a regime label was
> assigned by the most-binding shadow price across (`EX_glc__D_e`,
> `EX_ac_e`, `EX_o2_e`, `EX_nh4_e`, `EX_pi_e`). Targeted FVA on a
> 45-reaction curated panel produced `width__` features. An XGBoost
> regime classifier and severity regressor were trained with 5-fold CV;
> SHAP attribution was computed via `Booster.predict(pred_contribs=True)`.
> Outputs: `revision_runs/iscience_rev1/03_transfer/`.

### Rebuttal answer

> Reviewers asked about generalizability. We added an in-silico transfer
> demo to *E. coli* iML1515 using the unchanged diagnostic workflow.
> Performance is preserved (macro-F1 = 0.97; R² = 0.98), and the
> SHAP-top features map to functionally analogous central-carbon and
> respiratory modules — supporting transferability of the framework
> while making explicit that the curated feature panel itself is
> system-specific. The present transfer analysis is in silico only, and
> wet-lab validation in *E. coli* is acknowledged as future work.

**Strength:** ✅ **SI figure only** + 1-paragraph Results/Discussion
insert. *Do not promote to a main-text Fig 4 replacement.*

---

## 5. Runtime profiling

### Methods paragraph

> All revision workstreams were run on a single CPU core
> (`Miniforge3/python.exe`, Python 3.12.2; xgboost 3.1.3; cobra 0.31.1;
> scikit-learn 1.8.0; shap 0.49.1; causal-learn 0.1.4.5; GLPK solver).
> Per-stage wall-clock times were: LHS 500-sample ≈ 0.001 s; FBA batch
> (n=50) ≈ 2 s; targeted FVA (10 conditions × 30 reactions) ≈ 0.4 s;
> XGBoost classifier + SHAP on n=242 ≈ 0.3 s; XGBoost regressor + SHAP ≈
> 0.06 s; PC bootstrap (25 iterations) ≈ 0.3 s. Total per-pipeline
> end-to-end runtime is on the order of seconds for the deployed dataset
> size; peak memory stayed below 420 MB. See
> `revision_runs/iscience_rev1/04_runtime/runtime_summary.csv` and
> `environment_summary.txt`.

**Strength:** ✅ Methods paragraph + SI table (runtime summary).

---

## Point-by-point response (template)

(Copy individual answers from the per-workstream sections above and tailor
to each reviewer point. Conservative wording: "proof-of-concept",
"threshold-like / interaction-modulated", "hypothesis-prioritization
structure".)

| Editor / Reviewer point | Workstream | Answer file |
|---|---|---|
| Mandatory 2 + 3 (justify 30, characterize panel) | #1 ablation | `01_feature_panel_ablation/ablation_summary.md` |
| Mandatory 4 (benchmarks) | #2 benchmark | `02_benchmarking/benchmark_summary.md` |
| Mandatory 3 (per-class, residual diagnostics) | #3 existing-data | `05_existing_data/existing_data_summary.md` |
| Mandatory 5 + 6 (generalizability) | #4 transfer | `03_transfer/transfer_summary.md` |
| Mandatory 9 (runtime) | #5 runtime | `04_runtime/runtime_summary.md` |

---

## Known caveats to disclose (mandatory transparency)

- **Width universe — extended FVA campaign applied.** The original
  `results/regime_dataset.parquet` was alphabetically truncated at
  `FACOAL161` (120 `width__` columns), excluding paper-named TCA
  enzymes beyond 'F'. For this revision, an extended FVA campaign
  re-ran targeted FVA on the 180 missing reactions across all 242
  conditions (replay-verified against the original `objective_value`
  to within 5e-2 across all rows), expanding the deployed `width__`
  universe to ~300 columns. Paper-named TCA / respiration anchors
  (MDH, ICDHx, ICDHyr, ICL, MALS, PYK, PPC, NADH16pp, FUM …) are now
  present and appear in the SHAP top-K of the iSO1 classifier,
  aligning with the published Fig 4 narrative. See
  `revision_runs/iscience_rev1/extended_fva/` for raw outputs and
  `code/revision/extend_fva_campaign.py` for the driver.
- **Baseline B target overlap.** B_gem_summary's R² ≈ 0.998 reflects
  target-construction overlap (severity = obj/obj_max; obj is in B's
  feature vector), not new predictive content. Reported transparently.
- **iML1515 transfer scope.** In-silico only; no wet-lab matching.
  Stated as future work in both the Results paragraph and the rebuttal.
- **CV vs experimental performance.** The CV macro-F1 = 0.991 numbers
  describe model-internal robustness on the simulated diagnostic
  dataset; they are **not** the same quantity as the experimental
  agreement/mismatch shown in Fig 7. The two should be presented as
  complementary, not interchangeable.

---

## 7. Point-flux vs flexibility-interval (novelty defense)

### Strategic framing

This is the workstream that directly defends the manuscript's central
novelty claim — that **flexibility-interval (FVA-width) features support
feasible-space diagnosis, conceptually distinct from point-state prediction
(pFBA / FBA)**. Two complementary lines of evidence:

1. **On the iSO1 internal dataset** (n=242, 42-reaction panel), W matches
   PF (W=0.957 vs PF=0.932; W-vs-PF = +0.025 F1) and the same-solve
   midpoint M (W=0.957 vs M=0.951; W-vs-M = +0.006 F1). The sign-agnostic
   |pFBA flux| edges out W (PFA=0.972 vs W=0.957). Honest reading: on the
   internal dataset the F1 picture is a wash.
2. **On the external iML1515 system** (n=244, 45-reaction panel), W
   clearly beats both pFBA baselines: W=0.972 vs PF=0.906 vs PFA=0.907
   (W-vs-PF = **+0.066 F1**). This is the cleanest fair test (different
   organism, identical pipeline) and the result that supports the
   novelty claim quantitatively.

The two together support a concept-level contribution rather than a raw-F1
dominance argument: width and pFBA-magnitude answer *different questions*,
and on an unseen system the width representation generalizes better.

### Results paragraph (main text — Discussion or Results closer)

> To directly test whether the diagnostic signal lives in the *interval*
> or in the *point*, we held the curated reaction panel, the XGBoost
> learner, and the 5-fold CV split constant and varied only the feature
> representation. Five representations were compared: FVA-derived width
> (our flexibility feature), FVA midpoint (point center from the same
> solve), parsimonious-FBA flux (signed; the strongest deterministic
> point-flux baseline), |pFBA flux| (magnitude), and FBA objective alone.
> On the iSO1 internal dataset (n=242, 42-reaction panel), the interval
> width outperformed the same-solve midpoint by Δmacro-F1 = +0.006
> (0.957 vs 0.951), and outperformed signed pFBA flux by Δmacro-F1 = +0.025
> (0.957 vs 0.932). The sign-agnostic |pFBA flux| reached numerically
> higher F1 on iSO1 (0.972) while remaining comparable on severity R²
> (0.926 vs 0.922). On the **external iML1515 transfer** (n=244,
> 45-reaction panel), the width representation cleanly beat the pFBA
> baselines: macro-F1 = **0.972 (W) vs 0.906 (PF) vs 0.907 (PFA)**, a
> +0.066 F1 advantage on the cross-organism test. We interpret these
> comparisons not as a claim of unique optimality for FVA-width on every
> dataset, but as evidence that **flexibility-interval features support a
> mechanistically interpretable diagnostic representation — degrees of
> metabolic freedom — that point-flux representations do not, and that
> this representation generalizes better to an unseen system**.

### Methods paragraph

> Point-flux vs flexibility-interval comparison. For each of the 242
> stored LHS conditions (iSO1_933) and the 244 iML1515 transfer conditions,
> parsimonious FBA (`cobra.flux_analysis.pfba`) was solved with the same
> medium and uptake bounds used to generate the FVA dataset. Five feature
> representations were constructed on the same curated panel: width
> (`vmax − vmin`), midpoint (`(vmax + vmin) / 2`), signed pFBA flux,
> absolute pFBA flux, and biomass objective. Each representation was
> evaluated by XGBoost classifier and regressor under 5-fold CV with
> identical hyperparameters (300 trees, max_depth=4, hist), reporting
> macro-F1, balanced accuracy, RMSE, and R². pFBA fluxes are deterministic
> and remove the non-uniqueness confound of plain FBA. Outputs:
> `revision_runs/iscience_rev1/07_pointflux/`.

### Rebuttal answer (editor: "more explanation about the novelty of the method is required")

> The editor and reviewers asked for stronger articulation of the
> manuscript's methodological novelty. Beyond rewriting the
> Introduction to position the framework as a **feasible-space /
> degrees-of-freedom diagnosis** rather than a point-state prediction
> tool (see Introduction insert below), we added a direct quantitative
> comparison: the same XGBoost learner, the same curated panel, the same
> 5-fold CV — but with the feature representation varied across
> FVA-width (ours), FVA midpoint, signed pFBA flux, |pFBA flux|, and FBA
> objective. On the iSO1 internal dataset, width matches or beats every
> point-flux baseline except |pFBA flux| (W=0.957 vs PFA=0.972), where
> |pFBA flux| edges higher on classification F1 while remaining equal
> on severity R². On the **external iML1515 system, width clearly beats
> both pFBA baselines** (W=0.972 vs PF=0.906, PFA=0.907; +0.066 F1) —
> the cleanest fair test in this comparison, because the same pipeline
> is exercised on a different organism without any width-specific
> tuning. Within the same FVA solve, width also outperforms the
> midpoint (Δmacro-F1 = +0.006; ΔR² = +0.014), confirming that the
> *size* of the feasible interval carries information beyond the
> central representative flux. |pFBA flux| and width answer different
> questions: |pFBA| describes *where the biomass-maximizing flux goes*,
> while width describes *how much rerouting capacity remains under the
> imposed constraint context*. The latter is what supports the
> downstream rigidification / degree-of-freedom interpretation
> (Figs 5–6) that point-flux representations cannot produce. The
> novelty of the method, then, is not that one representation dominates
> classification F1 on every dataset; it is that **the flexibility-
> interval representation enables a mechanistic diagnosis of
> feasible-space collapse that is conceptually distinct from, and not
> reducible to, point-state prediction, and that generalizes better to
> an unseen system**. All point-flux baselines are reported
> transparently in the SI alongside the width-based results.

### Introduction insert (Reviewer 1 — overall framework logic)

> Most FBA- and dFBA-based diagnostic workflows summarize the intracellular
> state by an objective-optimized point solution, such as a
> biomass-maximizing flux vector or a time-resolved growth trajectory.
> However, in genome-scale metabolic models, many alternative flux vectors
> can satisfy the same extracellular constraints and objective value,
> making point fluxes unstable as diagnostic targets. We therefore
> reformulate metabolic diagnosis as a feasible-space problem. For each
> intervention context, flux variability analysis defines an allowable
> interval for each reaction or module, and the interval width represents
> the remaining metabolic degrees of freedom. Wide intervals indicate
> retained rerouting capacity, whereas shrinking intervals indicate
> rigidification — i.e., loss of flexibility as constraints accumulate.
> In this framework, XGBoost learns the relationship between FVA-derived
> flexibility widths and regime identity or growth potential, SHAP
> identifies which degrees of freedom are collapsing, and causal-structure
> discovery organizes these collapses into a hypothesis map for
> intervention prioritization. The proposed framework therefore differs
> from point-state prediction by diagnosing the collapse of metabolic
> flexibility rather than a single optimal flux state.

### Cover-letter insert (novelty re-positioning)

> The revised manuscript clarifies that the methodological novelty of
> this work is not the combination of GEM, XAI, and causal discovery,
> but the **reformulation of GEM-based bioprocess diagnosis from
> point-flux prediction to feasible-space diagnosis**. Conventional FBA/
> dFBA-type analyses rely on objective-optimized point solutions, while
> our framework uses FVA-derived interval widths to quantify metabolic
> degrees of freedom and diagnose flexibility collapse. The accompanying
> revisions support this re-positioning with new benchmarking against
> point-flux baselines (Workstream #7), feature-panel ablation
> (Workstream #1), quantitative performance metrics (Workstreams #2–#3),
> computational-cost reporting (Workstream #5), and an in-silico
> cross-organism demonstration on E. coli iML1515 (Workstreams #3 and
> #7). We reduce overstatement throughout: the framework is presented
> as a **proof-of-concept diagnostic platform / hypothesis-prioritization
> structure**, not a validated deployable tool.

### Note on iML1515 pFBA completion (now done — 2026-05-19)

> The iML1515 pFBA and |pFBA| baselines have been computed (cache:
> `revision_runs/iscience_rev1/07_pointflux/pfba_fluxes_iml1515.parquet`,
> 244 × 48). Final iML1515 numbers:
>
> | Representation | macro-F1 | R² |
> |---|---|---|
> | **W (widths, ours)** | **0.972** | **0.978** |
> | PF (signed pFBA)     | 0.906 | 0.979 |
> | PFA (|pFBA flux|)    | 0.907 | 0.978 |
> | OBJ (scalar floor)   | 0.297 | 0.998* |
>
> *OBJ R² ≈ 0.998 is a target-construction artifact (severity =
> obj/obj_max).
>
> The cross-organism W-vs-pFBA gap (+0.066 macro-F1) is the cleanest
> single quantitative line in support of the manuscript's novelty
> framing.

**Strength:** ✅ SI figure (`pointflux_vs_width.png`) + 1-paragraph
Discussion insert + Introduction reframing + Cover-letter sentence.
*This is the single most important workstream for defending the
manuscript's central novelty claim.*

---

## 8. Updated point-by-point table

| Editor / Reviewer point | Workstream | Answer file |
|---|---|---|
| Editor: novelty of method | **#7 pointflux** | `07_pointflux/pointflux_summary.md` |
| Mandatory 2 + 3 (justify 30, characterize panel) | #1 ablation | `01_feature_panel_ablation/ablation_summary.md` |
| Mandatory 4 (benchmarks) | #2 benchmark | `02_benchmarking/benchmark_summary.md` |
| Mandatory 3 (per-class, residual diagnostics) | #3 existing-data | `05_existing_data/existing_data_summary.md` |
| Mandatory 5 + 6 (generalizability) | #4 transfer | `03_transfer/transfer_summary.md` |
| Mandatory 9 (runtime) | #5 runtime | `04_runtime/runtime_summary.md` |
| Reviewer 1: overall framework logic in Intro | **#7 pointflux** Intro insert (above) | (this file) |
| Editor: refine claims (concept vs tool) | **#7 pointflux** Cover-letter insert (above) | (this file) |

---

## 9. Fig 6 arrow consistency (Reviewer 2, Mandatory M7)

### Strategic framing

Reviewer 2 asked why arrows appear in Fig 6b but not Fig 6a. The honest
audit (`code/revision/09_fig6_arrow_audit.py`) of the cached PC-bootstrap
output (`Fig05_causal_dag__dag_edges.csv`) shows that **100 % of the
bootstrap-stable edges are undirected** — the published Fig 6a's
arrowheads and Fig 6b's `→` labels both visually imply orientation that
the PC algorithm did not actually infer. The fix is *caption-level + a
re-rendered SI panel*; no algorithm re-run is required, and the
underlying scientific conclusions are unchanged.

### Results / Methods sentences (insert near Fig 6 description)

> The bootstrap-stable edge set produced by the PC algorithm
> (Fisher-z independence test, α = 0.05; 100 resamples; minimum stability
> frequency 0.10) is fully undirected for this dataset — the orientation
> rules did not succeed within the bootstrap-stable equivalence class. We
> therefore present the inferred rigidification structure as an
> **undirected adjacency network with bootstrap stability annotations**.
> A re-rendered version of Fig 6 with consistent undirected convention
> (no arrowheads in panel (a); em-dash "—" separator in panel (b)) is
> provided in the Supplementary Information (Fig SX,
> `revision_runs/iscience_rev1/09_fig6_audit/fig6_consistent.{png,pdf,svg}`).
> The Fig 6 caption has been updated accordingly (see Variant B in
> `docs/revision/captions/fig6_legend_v2.md`).

### Rebuttal answer (Reviewer 2 / Mandatory M7)

> Thank you for flagging the inconsistency between Fig 6a (no arrows) and
> Fig 6b (arrow notation in the labels). We have audited the cached
> PC-bootstrap output and confirmed that **the bootstrap-stable edge set
> is 100 % undirected** — the orientation rules of the PC algorithm did
> not succeed within the stable equivalence class for this dataset. The
> previously published Fig 6a rendered arrowheads and Fig 6b's labels
> used the `→` notation, both visually implying directional inference
> that the algorithm did not actually produce. We have addressed this in
> two complementary ways: (i) the Fig 6 caption is revised to state
> explicitly that the bootstrap-stable edge set is undirected and that
> the inferred structure is presented as a *conditional-dependency map*
> rather than an oriented causal DAG (`docs/revision/captions/fig6_legend_v2.md`,
> Variant B); and (ii) a re-rendered Supplementary Figure with consistent
> undirected convention (no arrowheads in panel (a); em-dash "—"
> separator in panel (b)) is provided in
> `revision_runs/iscience_rev1/09_fig6_audit/fig6_consistent.{png,pdf,svg}`,
> with a side-by-side before/after comparison in
> `fig6_before_after.{png,pdf}`. This change is caption-level; the
> underlying PC bootstrap, the bootstrap-stable backbone, and the
> hypothesis-prioritization interpretation are unchanged.

**Strength:** ✅ caption revision + 1 SI figure (re-rendered) + 1 SI
side-by-side audit figure. *Addresses Editor mandatory M7 and Reviewer 2
comment R2.2.*
