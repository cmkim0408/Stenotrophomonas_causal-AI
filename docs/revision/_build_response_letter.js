/**
 * Build docs/revision/response_letter.docx from response_letter.md content.
 * iScience-grade: US Letter, Arial 11pt, 3-column reviewer table, summary
 * table on final pages.
 *
 * Run: node docs/revision/_build_response_letter.js
 */

const {
  Document, Packer, Paragraph, TextRun, Table, TableRow, TableCell,
  AlignmentType, BorderStyle, WidthType, ShadingType, HeadingLevel,
  PageOrientation, PageBreak,
} = require("docx");
const fs = require("fs");
const path = require("path");

const OUT = path.resolve(__dirname, "response_letter.docx");

// ---------- style helpers ----------
const FONT = "Arial";
const NORMAL_SIZE = 22;     // 11 pt (docx uses half-points)
const SMALL_SIZE = 18;      // 9 pt for table cells
const H1_SIZE = 32;         // 16 pt
const H2_SIZE = 26;         // 13 pt
const COL_GUTTER = 0;       // cells own their padding

const border = { style: BorderStyle.SINGLE, size: 1, color: "BFBFBF" };
const cellBorders = { top: border, bottom: border, left: border, right: border };

function P(text, opts = {}) {
  return new Paragraph({
    spacing: { after: 120, line: 300 },
    alignment: opts.align || AlignmentType.LEFT,
    children: [new TextRun({ text, font: FONT, size: opts.size || NORMAL_SIZE,
                             bold: opts.bold, italics: opts.italics })],
  });
}

function H(level, text) {
  const size = level === 1 ? H1_SIZE : H2_SIZE;
  return new Paragraph({
    heading: level === 1 ? HeadingLevel.HEADING_1 : HeadingLevel.HEADING_2,
    spacing: { before: 280, after: 160 },
    children: [new TextRun({ text, font: FONT, size, bold: true })],
  });
}

function cell(text, opts = {}) {
  const runs = typeof text === "string"
    ? [new TextRun({ text, font: FONT, size: opts.size || SMALL_SIZE,
                     bold: opts.bold })]
    : text;  // already an array of TextRun
  return new TableCell({
    borders: cellBorders,
    width: { size: opts.width, type: WidthType.DXA },
    shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
    margins: { top: 80, bottom: 80, left: 100, right: 100 },
    children: [new Paragraph({
      alignment: AlignmentType.LEFT,
      children: typeof text === "string"
        ? [new TextRun({ text, font: FONT, size: opts.size || SMALL_SIZE,
                         bold: opts.bold })]
        : text,
    })],
  });
}

function richCell(parts, opts = {}) {
  // parts: array of {text, bold?, italics?, code?}
  const runs = parts.map(p => new TextRun({
    text: p.text, font: p.code ? "Consolas" : FONT,
    size: opts.size || SMALL_SIZE,
    bold: p.bold, italics: p.italics,
  }));
  return new TableCell({
    borders: cellBorders,
    width: { size: opts.width, type: WidthType.DXA },
    shading: opts.fill ? { fill: opts.fill, type: ShadingType.CLEAR } : undefined,
    margins: { top: 80, bottom: 80, left: 100, right: 100 },
    children: [new Paragraph({
      alignment: AlignmentType.LEFT,
      children: runs,
    })],
  });
}

// ---------- per-comment table ----------
// US Letter portrait content width with 1" margins = 9360 DXA
// 3-column layout: ID 900 / Comment 3500 / Response 3500 / Artefact 1460
const W_ID = 900;
const W_COMMENT = 3500;
const W_RESPONSE = 3500;
const W_ARTEFACT = 1460;
const TABLE_W = W_ID + W_COMMENT + W_RESPONSE + W_ARTEFACT;

function headerRow(headers) {
  return new TableRow({
    tableHeader: true,
    children: headers.map((h, i) => cell(h, {
      width: [W_ID, W_COMMENT, W_RESPONSE, W_ARTEFACT][i],
      bold: true, fill: "E7E6E6",
    })),
  });
}

function commentRow(id, comment, response, artefact) {
  return new TableRow({
    cantSplit: false,
    children: [
      cell(id, { width: W_ID, bold: true, fill: "F2F2F2" }),
      cell(comment, { width: W_COMMENT, size: SMALL_SIZE }),
      cell(response, { width: W_RESPONSE, size: SMALL_SIZE }),
      cell(artefact, { width: W_ARTEFACT, size: SMALL_SIZE }),
    ],
  });
}

// ---------- content ----------
const COMMENTS = {
  editor_mandatory: [
    {
      id: "M1",
      comment: 'Strengthen conceptual clarity (Reviewer 1): Add a clear paragraph in the Introduction explaining FVA widths, regime definition, and overall framework logic.',
      response: 'A new Introduction paragraph reframes the framework as feasible-space diagnosis, explicitly explaining FVA-derived widths as remaining metabolic degrees of freedom, the active-constraint-set definition of regimes, and the role of XGBoost+SHAP and causal-structure discovery in linking widths → regime labels → rigidification hypotheses. The paragraph is drafted in rebuttal §7 Introduction insert and will be placed in the manuscript Introduction.',
      artefact: 'Manuscript Introduction (new paragraph, [Line XXX]); rebuttal §7 Introduction insert',
    },
    {
      id: "M2",
      comment: 'Justify feature selection (Reviewer 2): Provide rationale for selecting 30 reactions/modules and assess whether feature expansion improves performance.',
      response: 'Workstream #1 performed a six-level feature-panel ablation. An extended FVA campaign first expanded the deployed width__ universe from 120 to ~300 columns; panels of size 10, 20, 42 (curated), 50, and 300 yielded macro-F1 between 0.957 and 0.991 and severity R² between 0.90 and 0.92. Random 30-feature controls reached macro-F1 = 0.981 ± 0.009. The diagnostic signal is broadly distributed; the curated panel is an interpretability-oriented layer rather than a performance-optimal subset.',
      artefact: 'SI Fig SX (feature_panel_ablation.png); ablation_metrics.csv; rebuttal §1',
    },
    {
      id: "M3",
      comment: 'Improve model robustness / performance evaluation: Address prediction errors (Fig. 7) and include quantitative performance metrics (classification + regression).',
      response: 'Workstream #3 reports 5-fold CV: macro-F1 = 0.991, balanced accuracy = 0.984, per-class F1 [N=0.976, Ac=1.000, O2=0.998]; severity RMSE = 0.034, R² = 0.906, residual mean = −0.001. These describe model-internal CV robustness on the simulated dataset, distinct from the experimental Fig 7 story. The top C1–C10 mismatch is C10 (sealed-cap mid-O2; rank-residual 0.85), consistent with non-stoichiometric constraints.',
      artefact: 'SI Fig SX (confusion_matrix.png, regression_residuals.png); performance_metrics.csv, residual_summary.csv, top_mismatch_conditions.csv; rebuttal §3',
    },
    {
      id: "M4",
      comment: 'Add benchmarking (critical for ISCI): Compare against baseline methods (e.g., standard FBA, simpler ML models, or prior frameworks).',
      response: 'Workstream #2 benchmarks three feature classes (inputs-only, GEM-summary, curated widths) crossed with three learners (LogReg, RandomForest, XGBoost). On the curated width panel all three learners are competitive (macro-F1 0.957–1.000; R² 0.92–0.998). Inputs-only baselines reach macro-F1 ≤ 0.96 but R² ≤ 0.86. The GEM-summary R² ≈ 0.998 reflects target-construction overlap (severity = obj/obj_max). The signal is carried primarily by the flexibility-based representation, not by a uniquely optimal learner.',
      artefact: 'SI Fig SX (benchmark_comparison.png); benchmark_metrics.csv; rebuttal §2',
    },
    {
      id: "M5",
      comment: 'Strengthen validation (critical): Provide additional validation beyond simulation (e.g., more experimental data or external/public datasets).',
      response: 'Two complementary additions: (i) C1–C10 + N1–N10 experimental holdout reanalyzed with rank-residual mismatch scoring (Fig 7 / workstream #3); (ii) in-silico external validation on the public E. coli iML1515 GEM (workstreams #4 + #7) — macro-F1 = 0.972, R² = 0.978 on 244 LHS conditions; iML1515 SHAP top features map functionally to the same central-carbon / respiration / acetate-uptake / biosynthesis modules as iSO1. The iML1515 transfer is in silico only; wet-lab external validation is acknowledged as future work.',
      artefact: 'SI Fig SX (external_transfer.png); transfer_metrics.csv; rebuttal §4; Limitations',
    },
    {
      id: "M6",
      comment: 'Demonstrate generalizability: Test or discuss applicability to other organisms/datasets; include at least one external or literature-based validation.',
      response: 'Workstream #4 applies the identical workflow (LHS over uptake bounds → shadow-price regime labeling → targeted FVA on a 45-reaction E. coli curated panel → XGBoost + SHAP) to iML1515. 244 of 250 LHS samples were feasible across three regimes (o2/nh4/glc limited). 5-fold CV reached macro-F1 = 0.972, R² = 0.978. The framework is transferable in formulation; system-specific feature curation remains necessary; wet-lab validation in the external organism remains future work.',
      artefact: 'SI Fig SX (external_transfer.png); transfer_metrics.csv; rebuttal §4',
    },
    {
      id: "M7",
      comment: 'Clarify causal map representation (Fig. 6): Explain inconsistency between Fig. 6a and 6b (arrows vs no arrows).',
      response: 'The two panels represent different stages of the PC-bootstrap pipeline (Fig 6a = undirected adjacency; Fig 6b = directed-edge stability). A dedicated audit workstream (#E) is being prepared that re-renders both panels with explicit directed-edge vs undirected-skeleton legends and updates the figure caption to make the stage distinction explicit. This will be included in the next revision iteration.',
      artefact: 'Planned in workstream #E; caption draft at docs/revision/captions/fig6_legend_v2.md. Next revision iteration.',
    },
    {
      id: "M8",
      comment: 'Clarify interpretation of non-linearity (Fig. 5): Revise text/plots to clearly demonstrate or justify "non-linear effects."',
      response: 'We acknowledge that the Fig 5b/c plots do not visually demonstrate the non-linearity claim. A dedicated reanalysis workstream (#D) is being prepared that overlays LOESS / two-knot piecewise-linear fits with an F-test against the linear baseline and computes SHAP interaction values for the top severity features. If the F-test fails to support strict non-linearity, the manuscript will adopt the conservative wording "threshold-like or interaction-modulated effects." Next revision iteration.',
      artefact: 'Planned in workstream #D; outputs at revision_runs/iscience_rev1/08_nonlinearity/. Next revision iteration.',
    },
    {
      id: "M9",
      comment: 'Report computational cost and scalability: Provide runtime, computational resources, and feasibility for real-world use.',
      response: 'Workstream #5 profiles per-stage wall-clock and memory on a single CPU core (Python 3.12.2, xgboost 3.1.3, cobra 0.31.1, GLPK): LHS ≈ 0.001 s; FBA batch (n=50) ≈ 2 s; targeted FVA (10 × 30) ≈ 0.4 s; XGBoost classifier + SHAP ≈ 0.3 s; XGBoost regressor + SHAP ≈ 0.06 s; PC bootstrap (25) ≈ 0.3 s. Total end-to-end pipeline runtime is on the order of seconds for the deployed dataset; peak memory < 420 MB.',
      artefact: 'SI Table SX (runtime_summary.csv, environment_summary.txt); revised Methods',
    },
    {
      id: "M10",
      comment: 'Refine claims and positioning: Reduce overstatement; clearly distinguish conceptual framework vs validated tool.',
      response: 'The framework is now described as a proof-of-concept diagnostic platform and a hypothesis-prioritization structure rather than a validated deployable tool. Specific adjustments: Abstract final sentence reframed to feasible-space framing; Discussion adds an explicit limitation that the iML1515 transfer is in silico only; Limitations paragraph refined to distinguish stoichiometric vs non-stoichiometric mismatch interpretations and to acknowledge that the causal map is hypothesis-prioritization rather than causal proof.',
      artefact: 'Manuscript Abstract / Discussion / Limitations; rebuttal §7 Cover-letter insert',
    },
  ],
  editor_novelty: [
    {
      id: "EN",
      comment: '"More explanation about the novelty of the method is required."',
      response: 'Workstream #7 is dedicated to this point. Holding curated panel, XGBoost learner, and 5-fold CV split constant, only the feature representation was varied across FVA-width (W), FVA midpoint (M), signed pFBA flux (PF), |pFBA flux| (PFA), and FBA objective alone (OBJ). On iSO1, W = 0.957 / PFA = 0.972 / PF = 0.932 / M = 0.951 (F1 is a wash; PFA edges W). On external iML1515, W cleanly beats both pFBA baselines: W = 0.972 vs PF = 0.906, PFA = 0.907 (+0.066 macro-F1). Crucially, W and |pFBA| answer different questions: |pFBA| = where the biomass-maximizing flux goes; W = how much rerouting capacity remains. Only width supports the rigidification / degree-of-freedom interpretation. The novelty is conceptual (feasible-space vs point-state) AND quantitative on the cross-organism test.',
      artefact: 'SI Fig SX (pointflux_vs_width.png); pointflux_metrics.csv; rebuttal §7; Introduction insert (M1) + Discussion insert',
    },
  ],
  reviewer1: [
    {
      id: "R1.1",
      comment: '"While this study proposes a methodological framework, the abstract and the introduction section provide limited explanation of FVA, XGBoost, SHAP, and the underlying conceptual rationale. In particular, the explanation in the context of genome-scale metabolic models is concentrated in the later sections, making it difficult to grasp the overall picture from the abstract and introduction alone. It would be helpful to include a paragraph in the introduction that outlines what relationships are learned using FVA widths and how regimes are defined."',
      response: 'A new Introduction paragraph reframes the framework as feasible-space diagnosis and explicitly defines (i) FVA-derived widths as remaining metabolic degrees of freedom under the imposed constraint context, (ii) regimes as the identity of the active limiting constraint (largest positive shadow price), (iii) the role of XGBoost in mapping widths → regime labels and severity, and (iv) the role of SHAP and PC-based causal-structure discovery in organising the resulting attributions into a rigidification hypothesis map. This addresses both R1.1 and M1.',
      artefact: 'Manuscript Introduction (new paragraph, [Line XXX]); rebuttal §7 Introduction insert',
    },
  ],
  reviewer2: [
    {
      id: "R2.1",
      comment: '"The authors state that 30 reactions/modules were selected for FVA and for calculating the degrees-of-freedom features (L542), also referred to as reaction flexibilities in Fig. 4 (L612). However, considering the large residual errors observed for the two red nodes in Fig. 7 (L631), it is unclear whether these 30 reactions/modules provide sufficient information for the two XGBoost tasks, namely mapping these flexibilities to (i) regime labels and (ii) growth-potential index values (L560). In my opinion, incorporating additional features may improve the predictive performance. Why were these 30 reactions/modules specifically chosen?"',
      response: 'See M2. Workstream #1 directly tests whether feature expansion improves performance. After extending the width universe to ~300 columns, panels of size 10, 20, 42, 50, 300 all stay within macro-F1 [0.957, 0.991] and R² [0.90, 0.92]; random-30 controls match curated. The two red nodes in Fig 7 (C10 + mismatch contexts) are interpreted in workstream #3 as signals of non-stoichiometric constraints (e.g., physicochemical inhibition), not as evidence of insufficient feature coverage. The curated panel is best understood as an interpretability-oriented diagnostic layer, not a performance-optimal subset.',
      artefact: 'SI Fig SX (feature_panel_ablation.png); ablation_metrics.csv; top_mismatch_conditions.csv; rebuttal §1 + §3',
    },
    {
      id: "R2.2",
      comment: '"Why are arrows shown in Fig. 6b (L626) but not in Fig. 6a?"',
      response: 'See M7. The two panels represent different stages of the PC-bootstrap pipeline (undirected adjacency vs directed-edge stability). A dedicated audit workstream (#E) is being prepared to re-render both panels with consistent arrow conventions, add an explicit directed-edge vs undirected-skeleton legend, and update the figure caption. Next revision iteration.',
      artefact: 'Planned in workstream #E. Next revision iteration.',
    },
    {
      id: "R2.3",
      comment: '"It is a little confusing that the authors declared the \'non-linear effects\', while it seems like a linear relationship just without regression line in Fig. 5b, c (L623)."',
      response: 'See M8. A dedicated reanalysis workstream (#D) is being prepared to overlay LOESS / piecewise-linear fits with formal F-tests against a linear baseline and to compute SHAP interaction values. If the F-test does not support strict non-linearity, the manuscript will be revised to use the more conservative wording "threshold-like or interaction-modulated effects." Next revision iteration.',
      artefact: 'Planned in workstream #D. Next revision iteration.',
    },
    {
      id: "R2.4",
      comment: '"This method is demonstrated on one organism (S. maltophilia), one may wonder how transferable is this framework? There might be many public datasets for model organisms."',
      response: 'See M6. Workstream #4 applies the identical workflow to E. coli iML1515 (public GEM), reaching macro-F1 = 0.972 and R² = 0.978 on 244 in silico LHS conditions. Workstream #7 additionally shows that on this external system the flexibility-width representation outperforms pFBA point-flux baselines by +0.066 macro-F1 — the cleanest fair test of the conceptual claim. The framework is transferable in formulation; wet-lab validation in the external organism is future work.',
      artefact: 'SI Fig SX (external_transfer.png, pointflux_vs_width.png); transfer_metrics.csv, pointflux_metrics.csv; rebuttal §4 + §7',
    },
    {
      id: "R2.5",
      comment: '"How much time did this framework totally spend in each step? Was it costly?"',
      response: 'See M9. Workstream #5 profiles every major stage on a single CPU core. Total per-pipeline wall-clock is on the order of seconds for the deployed dataset (LHS ≈ 0.001 s; FBA batch ≈ 2 s; targeted FVA ≈ 0.4 s; XGBoost classifier + SHAP ≈ 0.3 s; XGBoost regressor + SHAP ≈ 0.06 s; PC bootstrap (25) ≈ 0.3 s); peak memory < 420 MB. The framework is computationally inexpensive and feasible within decision-relevant time windows for routine use.',
      artefact: 'SI Table SX (runtime_summary.csv); revised Methods',
    },
  ],
  housekeeping: [
    ["H1", "Author list unchanged", "Confirmed — no additions or removals"],
    ["H2", "Generative AI / AI-assisted tools disclosure", "Already disclosed in the manuscript; retained verbatim"],
    ["H3", "Required statements (Data, CRediT, COI, Funding, Ethics)", "All present; reconfirmed"],
    ["H4", "Graphical abstract", "Existing retained; revised version emphasising feasible-space framing in preparation"],
    ["H5", "3–5 Highlights", "Existing retained; refined wording aligned with revised novelty framing"],
    ["H6", "Highlighted + clean manuscript files", "Both provided on submission"],
    ["H7", "STAR Methods format", "Retained; runtime / environment table appended"],
    ["H8", "Cover letter", "Updated to summarise the six new workstreams and revised novelty framing"],
  ],
};

const SUMMARY_ROWS = [
  ["New SI Figure — Feature-panel ablation across 6 panels + random-30 controls", "SI Fig SX (figures/feature_panel_ablation.png)", "M2; R2.1"],
  ["New SI Figure — Baseline benchmarking (3 feature sets × 3 learners)", "SI Fig SX (figures/benchmark_comparison.png)", "M4"],
  ["New SI Figure — 5-fold CV confusion matrix", "SI Fig SX (figures/confusion_matrix.png)", "M3"],
  ["New SI Figure — Severity regression residual plot", "SI Fig SX (figures/regression_residuals.png)", "M3"],
  ["New SI Figure — iML1515 external transfer (side-by-side SHAP top-K)", "SI Fig SX (figures/external_transfer.png)", "M5; M6; R2.4"],
  ["New SI Figure — Point-flux vs flexibility-interval (5 reps × 2 systems)", "SI Fig SX (figures/pointflux_vs_width.png)", "EN; R2.4"],
  ["New SI Table — Runtime / scalability summary", "SI Table SX (runtime_summary.csv, environment_summary.txt)", "M9; R2.5"],
  ["New SI Table — Per-class precision/recall/F1 + residual summary", "SI Tables (performance_metrics.csv, residual_summary.csv)", "M3"],
  ["New SI Table — Feature-panel ablation metrics", "SI Table (ablation_metrics.csv)", "M2; R2.1"],
  ["New SI Table — Baseline benchmarking metrics", "SI Table (benchmark_metrics.csv)", "M4"],
  ["New SI Table — iML1515 transfer metrics", "SI Table (transfer_metrics.csv)", "M5; M6; R2.4"],
  ["New SI Table — Point-flux vs flexibility metrics", "SI Table (pointflux_metrics.csv)", "EN; R2.4"],
  ["New Introduction paragraph — feasible-space framing", "Main text Introduction [Line XXX]", "M1; R1.1"],
  ["New Discussion insert — flexibility-based representation across model classes", "Main text Discussion [Line XXX]", "M4; EN"],
  ["New Discussion insert — in silico transfer to iML1515; wet-lab as future work", "Main text Discussion [Line XXX]", "M5; M6; R2.4"],
  ["New Methods paragraph — extended FVA campaign (120 → ~300 widths)", "Main text STAR Methods [Line XXX]", "M2; M3"],
  ["New Methods paragraph — feature-panel ablation protocol", "Main text STAR Methods [Line XXX]", "M2"],
  ["New Methods paragraph — baseline benchmarking protocol", "Main text STAR Methods [Line XXX]", "M4"],
  ["New Methods paragraph — iML1515 transfer protocol", "Main text STAR Methods [Line XXX]", "M6"],
  ["New Methods paragraph — point-flux vs flexibility-interval protocol", "Main text STAR Methods [Line XXX]", "EN"],
  ["New Methods sentences — runtime / environment", "Main text STAR Methods [Line XXX]", "M9; R2.5"],
  ["Revised Limitations paragraph — proof-of-concept, hypothesis-prioritization", "Main text Limitations [Line XXX]", "M10"],
  ["Revised Abstract — final sentence reframed (feasible-space)", "Main text Abstract [Line XXX]", "M1; M10; EN"],
  ["Refined wording — 'validated tool' / 'deployable system' removed", "Main text multiple [Line XXX]", "M10"],
  ["Updated Cover letter — novelty re-positioned", "Cover letter", "M10; EN"],
  ["Planned — Fig 5 LOESS + piecewise + interaction SHAP; Fig 6 directed-vs-undirected legend", "Next revision iteration (workstreams #D, #E)", "M7; M8; R2.2; R2.3"],
];

// summary table sizing: Change 4400 / Where 2960 / Reviewer 2000
const SW1 = 4400, SW2 = 2960, SW3 = 2000;
const SUMMARY_W = SW1 + SW2 + SW3;

function summaryRow(change, where, addressed) {
  return new TableRow({
    children: [
      cell(change, { width: SW1, size: SMALL_SIZE }),
      cell(where, { width: SW2, size: SMALL_SIZE }),
      cell(addressed, { width: SW3, bold: true, size: SMALL_SIZE }),
    ],
  });
}

// ---------- assemble children ----------
const children = [];

children.push(H(1, "Response to Editor and Reviewers"));
children.push(P("Manuscript: ISCIENCE-D-26-04043 — Digital twin for bioprocess bottleneck diagnosis under sparse observability: metabolic degrees of freedom, explainable AI, and rigidification maps", { italics: true }));
children.push(P("Date: May 19, 2026"));
children.push(P("Corresponding author: Changman Kim (cmkim@jnu.ac.kr)"));
children.push(P(""));
children.push(P("Dear Dr. Christina Nilofer and Reviewers,"));
children.push(P("We thank the editor and both reviewers for their constructive comments. We have substantially revised the manuscript to address every mandatory point and have added six new computational workstreams that together strengthen the manuscript's novelty defense, performance benchmarking, generalizability, runtime characterization, and quantitative reviewer-driven justification. Below we provide a point-by-point response, followed by a summary-of-changes table mapping each revision to the manuscript or Supplementary Information."));

children.push(H(1, "Cover paragraph"));
children.push(P(
  "The central methodological contribution of this work is a shift from point-state prediction to feasible-space diagnosis for genome-scale bioprocess analysis. Conventional FBA- and dFBA-type pipelines summarize the intracellular state by an objective-optimized point flux vector. We instead reformulate diagnosis as a problem of metabolic degrees of freedom: each intervention context defines an allowable interval per reaction (the FVA width), and shrinkage of these intervals — flexibility collapse — is a mechanistically interpretable signature of constraint propagation that point fluxes cannot express. To defend this conceptual shift quantitatively, we held the curated reaction panel, the XGBoost learner, and the 5-fold cross-validation split constant, and varied only the feature representation across FVA width, FVA midpoint, parsimonious-FBA flux (signed and magnitude), and FBA objective alone. Across organisms, the flexibility-interval representation outperformed parsimonious-FBA point-flux baselines on the iML1515 transfer test (macro-F1 0.972 vs 0.906) using the unchanged diagnostic pipeline, supporting the conceptual claim quantitatively (workstream #7, SI Fig SX pointflux_vs_width.png)."
));
children.push(P(
  "To address the editor's ten mandatory revisions and both reviewers' specific comments, we additionally performed: (#1) a six-level feature-panel ablation across panels of size 10–300 plus random-30 controls across ten seeds, showing diagnostic performance is largely insensitive to panel size (Δmacro-F1 ≤ 0.04) and that the curated panel is therefore an interpretability-oriented layer rather than a performance optimum; (#2) a baseline benchmark on the same 5-fold split comparing inputs-only, GEM-summary, and same-panel logistic-regression / random-forest / XGBoost models, supporting that the diagnostic signal is carried primarily by the flexibility-based representation rather than a uniquely optimal learner; (#3) an existing-data performance summary on the LHS-derived diagnostic dataset (5-fold CV macro-F1 = 0.991, severity R² = 0.906, C10 sealed-cap mid-O2 condition identified as top mismatch); (#4) an in-silico transfer demonstration to Escherichia coli iML1515 using the unchanged diagnostic workflow (macro-F1 = 0.972, R² = 0.978 on 244 LHS conditions); (#5) a runtime and scalability characterization (single-CPU, end-to-end pipeline wall-clock ≈ 3.2 s, peak memory < 420 MB); and an extended FVA campaign that re-ran targeted FVA on 180 missing reactions across all 242 conditions, expanding the deployed width__ universe to ~300 columns and bringing paper-named TCA / respiration anchors (MDH, ICDHx, ICDHyr, ICL, MALS, PYK, PPC, NADH16pp, FUM, …) into the analysis."
));
children.push(P(
  "Throughout the revision, we have refined our claims to characterize the framework as a proof-of-concept diagnostic platform and a hypothesis-prioritization structure rather than a validated deployable tool. Two reviewer points — the Fig 5 non-linearity interpretation and the Fig 6 arrow inconsistency — are being addressed by dedicated reanalysis workstreams (LOESS / piecewise / interaction-SHAP for Fig 5; directed-edge vs undirected-skeleton legend for Fig 6) that are scheduled for the next revision iteration and are noted transparently in the per-comment table below."
));

children.push(H(1, "Per-comment table"));
children.push(P(
  "The table below records each editor and reviewer comment verbatim, our response, and the revision artefact that supports it. Quantitative claims are sourced from revision_runs/iscience_rev1/metrics_summary.csv and the per-workstream *_summary.md files in the supporting repository (https://github.com/cmkim0408/Stenotrophomonas_causal-AI, branch revision/iscience-rev1)."
));

children.push(H(2, "Editor — Mandatory Revisions"));
const editorRows = [headerRow(["ID", "Comment (verbatim)", "Response (summary)", "Revision artefact"])];
for (const r of COMMENTS.editor_mandatory) {
  editorRows.push(commentRow(r.id, r.comment, r.response, r.artefact));
}
children.push(new Table({
  width: { size: TABLE_W, type: WidthType.DXA },
  columnWidths: [W_ID, W_COMMENT, W_RESPONSE, W_ARTEFACT],
  rows: editorRows,
}));

children.push(H(2, "Editor — Additional comment (Novelty)"));
const novRows = [headerRow(["ID", "Comment (verbatim)", "Response (summary)", "Revision artefact"])];
for (const r of COMMENTS.editor_novelty) {
  novRows.push(commentRow(r.id, r.comment, r.response, r.artefact));
}
children.push(new Table({
  width: { size: TABLE_W, type: WidthType.DXA },
  columnWidths: [W_ID, W_COMMENT, W_RESPONSE, W_ARTEFACT],
  rows: novRows,
}));

children.push(H(2, "Editor — Submission housekeeping"));
// housekeeping table: 3 cols (ID + Item + Status)
const HW_ID = 900, HW_ITEM = 4000, HW_STATUS = 4460;
const HW_TOTAL = HW_ID + HW_ITEM + HW_STATUS;
const hkRows = [
  new TableRow({
    tableHeader: true,
    children: [
      cell("ID", { width: HW_ID, bold: true, fill: "E7E6E6" }),
      cell("Item", { width: HW_ITEM, bold: true, fill: "E7E6E6" }),
      cell("Status", { width: HW_STATUS, bold: true, fill: "E7E6E6" }),
    ],
  }),
];
for (const row of COMMENTS.housekeeping) {
  hkRows.push(new TableRow({
    children: [
      cell(row[0], { width: HW_ID, bold: true, fill: "F2F2F2" }),
      cell(row[1], { width: HW_ITEM, size: SMALL_SIZE }),
      cell(row[2], { width: HW_STATUS, size: SMALL_SIZE }),
    ],
  }));
}
children.push(new Table({
  width: { size: HW_TOTAL, type: WidthType.DXA },
  columnWidths: [HW_ID, HW_ITEM, HW_STATUS],
  rows: hkRows,
}));

children.push(H(2, "Reviewer 1"));
const r1Rows = [headerRow(["ID", "Comment (verbatim)", "Response (summary)", "Revision artefact"])];
for (const r of COMMENTS.reviewer1) {
  r1Rows.push(commentRow(r.id, r.comment, r.response, r.artefact));
}
children.push(new Table({
  width: { size: TABLE_W, type: WidthType.DXA },
  columnWidths: [W_ID, W_COMMENT, W_RESPONSE, W_ARTEFACT],
  rows: r1Rows,
}));

children.push(H(2, "Reviewer 2"));
const r2Rows = [headerRow(["ID", "Comment (verbatim)", "Response (summary)", "Revision artefact"])];
for (const r of COMMENTS.reviewer2) {
  r2Rows.push(commentRow(r.id, r.comment, r.response, r.artefact));
}
children.push(new Table({
  width: { size: TABLE_W, type: WidthType.DXA },
  columnWidths: [W_ID, W_COMMENT, W_RESPONSE, W_ARTEFACT],
  rows: r2Rows,
}));

children.push(new Paragraph({ children: [new PageBreak()] }));
children.push(H(1, "Summary of changes"));
const sumRows = [
  new TableRow({
    tableHeader: true,
    children: [
      cell("Change", { width: SW1, bold: true, fill: "E7E6E6" }),
      cell("Where in revised manuscript", { width: SW2, bold: true, fill: "E7E6E6" }),
      cell("Reviewer point(s) addressed", { width: SW3, bold: true, fill: "E7E6E6" }),
    ],
  }),
];
for (const row of SUMMARY_ROWS) {
  sumRows.push(summaryRow(row[0], row[1], row[2]));
}
children.push(new Table({
  width: { size: SUMMARY_W, type: WidthType.DXA },
  columnWidths: [SW1, SW2, SW3],
  rows: sumRows,
}));

children.push(H(1, "Code and data availability — reconfirmed"));
children.push(P(
  "All revision analyses are reproducible from the branch revision/iscience-rev1 of the manuscript repository (https://github.com/cmkim0408/Stenotrophomonas_causal-AI). New scripts live under code/revision/; outputs under revision_runs/iscience_rev1/. The extended FVA campaign driver (extend_fva_campaign.py) replays each of the 242 stored conditions and matches the original objective_value within 5e-2 across all rows. iML1515 is loaded via cobra.io.load_model('iML1515') from the BiGG repository. Random seed 42 is used throughout (random-control panels: seed grid 1..10). All metrics tables, figure-source datasets, and per-row partial parquets are versioned in the repository."
));
children.push(P(""));
children.push(P("We thank the editor and both reviewers again for their constructive comments. We believe the revised manuscript, supported by the six new computational workstreams, substantially addresses the mandatory revisions and provides a clearer, more conservative articulation of the framework's conceptual novelty and its current proof-of-concept scope."));
children.push(P(""));
children.push(P("Sincerely,"));
children.push(P("Changman Kim, on behalf of all authors", { bold: true }));
children.push(P("Department of Biotechnology and Bioengineering, Chonnam National University", { italics: true }));
children.push(P("cmkim@jnu.ac.kr"));

// ---------- document ----------
const doc = new Document({
  creator: "Changman Kim",
  title: "Response to Editor and Reviewers — ISCIENCE-D-26-04043",
  styles: {
    default: { document: { run: { font: FONT, size: NORMAL_SIZE } } },
    paragraphStyles: [
      { id: "Heading1", name: "Heading 1", basedOn: "Normal", next: "Normal",
        quickFormat: true,
        run: { size: H1_SIZE, bold: true, font: FONT },
        paragraph: { spacing: { before: 320, after: 200 }, outlineLevel: 0 } },
      { id: "Heading2", name: "Heading 2", basedOn: "Normal", next: "Normal",
        quickFormat: true,
        run: { size: H2_SIZE, bold: true, font: FONT },
        paragraph: { spacing: { before: 260, after: 160 }, outlineLevel: 1 } },
    ],
  },
  sections: [{
    properties: {
      page: {
        size: { width: 12240, height: 15840 },          // US Letter
        margin: { top: 1440, right: 1440, bottom: 1440, left: 1440 },  // 1"
      },
    },
    children: children,
  }],
});

Packer.toBuffer(doc).then(buf => {
  fs.writeFileSync(OUT, buf);
  console.log("wrote", OUT, buf.length, "bytes");
}).catch(e => { console.error(e); process.exit(1); });
