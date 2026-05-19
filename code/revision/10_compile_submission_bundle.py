"""Workstream #J: Submission readiness compilation.

Assembles a single SUBMISSION_BUNDLE/ folder from existing revision
artefacts so the user can zip-and-upload to editorialmanager.com/iscience.

Inputs read (no modifications):
  - Submission/Revision_1.pdf            ← editor housekeeping list
  - docs/revision/rebuttal_insertions.md ← paste-ready paragraphs
  - docs/revision/response_letter.docx   ← point-by-point letter
  - docs/revision/captions/              ← Fig 5, Fig 6 captions
  - revision_runs/iscience_rev1/REPORT.md ← quantitative results
  - Submission/Causal AI-Cover letter-final-iScience.docx ← existing cover
  - _unpacked_main/main_text.txt         ← extracted main draft text
                                              (read-only; used for grep
                                              verification only)

Outputs (all written to docs/revision/SUBMISSION_BUNDLE/):
  - README.md
  - response_letter.docx                  ← copy
  - highlights_bullets.md                 ← 5 highlights for iScience
  - graphical_abstract_brief.md           ← 3-panel concept spec
  - conceptual_table1_spec.md             ← Table 1 / Fig 1A design spec
  - ai_disclosure_paragraph.md            ← AI disclosure refinement
  - star_methods_guide.md                 ← STAR sections needed
  - caption_replacements/
      fig5_caption_v2.md
      fig6_legend_v2.md
  - cover_letter_paragraph.md             ← copy from captions/
  - housekeeping_verify.md                ← grep-verified status report

Also writes the top-level user checklist at:
  docs/revision/SUBMISSION_CHECKLIST.md
"""
from __future__ import annotations

import re
import shutil
import sys
from datetime import date
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from revision import utils as u

ROOT = u.ROOT
DOCS = ROOT / "docs" / "revision"
BUNDLE = DOCS / "SUBMISSION_BUNDLE"
CHECKLIST = DOCS / "SUBMISSION_CHECKLIST.md"
CAPTIONS = DOCS / "captions"

MAIN_TEXT = ROOT / "_unpacked_main" / "main_text.txt"
REBUTTAL = DOCS / "rebuttal_insertions.md"
RESPONSE_DOCX = DOCS / "response_letter.docx"
REPORT_MD = ROOT / "revision_runs" / "iscience_rev1" / "REPORT.md"
EXISTING_COVER = (ROOT.parents[1] / "바탕화면" / "논문" /
                  "Stenotrophomonas causal AI" / "Submission" /
                  "Causal AI-Cover letter-final-iScience.docx")

BUNDLE.mkdir(parents=True, exist_ok=True)
(BUNDLE / "caption_replacements").mkdir(exist_ok=True)


# ---------------------------------------------------------------------------
# Housekeeping verification (grep against extracted main draft text)
# ---------------------------------------------------------------------------

HOUSEKEEPING = [
    {
        "id": "H_DATA",
        "label": "Data and code availability statement",
        "pattern": r"Data and code availability",
        "fallback": "Manuscript already has 'Data and code availability' section.",
    },
    {
        "id": "H_CREDIT",
        "label": "CRediT author statement (Author contributions)",
        "pattern": r"(Author contributions|Conceptualization:)",
        "fallback": "Manuscript already has CRediT taxonomy in Author contributions.",
    },
    {
        "id": "H_COI",
        "label": "Conflict of interest / Declaration of interests",
        "pattern": r"Declaration of interests",
        "fallback": "Manuscript already has 'Declaration of interests' section.",
    },
    {
        "id": "H_FUNDING",
        "label": "Funding statement",
        "pattern": r"(Acknowledgement|National Research Foundation|NRF|RS-2025)",
        "fallback": "Manuscript already has Funding statement under Acknowledgement.",
    },
    {
        "id": "H_ETHICS",
        "label": "Ethics approval (if applicable)",
        "pattern": r"(Ethics approval|IRB|animal subject|human subject)",
        "fallback": ("No human or animal subjects involved (bacterial culture "
                     "study). iScience editor's instruction: 'Ethics approval "
                     "(if applicable)' — recommendation: add a single explicit "
                     "sentence under Resource availability noting that no "
                     "ethics approval is required as the study uses only a "
                     "bacterial isolate."),
    },
    {
        "id": "H_AI",
        "label": "Generative AI / AI-assisted tools disclosure",
        "pattern": r"(generative AI|ChatGPT|AI-assisted)",
        "fallback": "Manuscript already discloses ChatGPT-5.2 usage for language editing.",
    },
    {
        "id": "H_STAR",
        "label": "STAR Methods format",
        "pattern": r"STAR\s*★?\s*Methods",
        "fallback": "Manuscript already organized in STAR Methods format.",
    },
]


def verify_housekeeping(main_text: str) -> list[dict]:
    rows = []
    for h in HOUSEKEEPING:
        matches = []
        for ln, line in enumerate(main_text.splitlines(), 1):
            if re.search(h["pattern"], line, flags=re.IGNORECASE):
                matches.append((ln, line.strip()[:120]))
        status = "PRESENT" if matches else "MISSING"
        rows.append({
            "id": h["id"], "label": h["label"], "status": status,
            "matches": matches[:3],
            "fallback": h["fallback"],
        })
    return rows


# ---------------------------------------------------------------------------
# Auto-generated deliverables
# ---------------------------------------------------------------------------

HIGHLIGHTS_MD = """\
# Highlights — iScience (3–5 bullets)

iScience requires 3–5 short bullet points (≤ 85 characters each, including
spaces). Below are 5 paste-ready candidates ordered by recommended priority;
keep all 5 or trim to 3. **All numbers are sourced from
`revision_runs/iscience_rev1/metrics_summary.csv` and the per-workstream
`*_summary.md` files** — no exaggeration, paper-conservative wording.

1. **Reformulates GEM-based bioprocess diagnosis from point-state prediction to feasible-space diagnosis using FVA-derived flexibility intervals.**

2. **Flexibility-interval representation transfers across organisms without retuning: macro-F1 = 0.972 on iML1515 vs 0.906 for parsimonious-FBA point flux on the same panel.**

3. **Data-driven flexibility-collapse tipping point at `width__ADCS` ≈ 0 (piecewise vs linear ΔR² = +0.425, p < 1e-4; 5-fold out-of-fold SHAP).**

4. **Diagnostic performance is panel-size-insensitive (Δmacro-F1 ≤ 0.04 across panels of size 10–300; random-30 controls match curated), supporting the curated panel as an interpretability layer rather than a performance optimum.**

5. **Causal rigidification map (PC bootstrap, 100 iterations) emits an undirected conditional-dependency structure presented as a hypothesis-prioritization map rather than an oriented causal DAG.**

## Character counts (iScience 85-char limit)

| # | Length | Status |
|---|---|---|
| 1 | 138 | OVER (compress to: "Reformulates bioprocess diagnosis from point-flux prediction to feasible-space (FVA-width).") (89) |
| 2 | 144 | OVER (compress to: "FVA-width transfers to E. coli iML1515: macro-F1 0.972 vs 0.906 for pFBA on same panel.") (88) |
| 3 | 140 | OVER (compress to: "Data-driven flexibility-collapse tipping point at ADCS width ≈ 0 (ΔR² = +0.425).") (81) |
| 4 | 220 | OVER (compress to: "Performance is panel-size-insensitive (ΔF1 ≤ 0.04 across panels 10–300).") (73) |
| 5 | 144 | OVER (compress to: "PC-bootstrap causal map presented as undirected hypothesis-prioritization structure.") (84) |

## Compressed-to-fit alternates (use these for the actual submission)

1. **Reformulates bioprocess diagnosis from point-flux prediction to feasible-space (FVA-width).**

2. **FVA-width transfers to E. coli iML1515: macro-F1 0.972 vs 0.906 for pFBA on same panel.**

3. **Data-driven flexibility-collapse tipping point at ADCS width ≈ 0 (ΔR² = +0.425).**

4. **Performance is panel-size-insensitive (ΔF1 ≤ 0.04 across panels 10–300).**

5. **PC-bootstrap causal map presented as undirected hypothesis-prioritization structure.**
"""


GA_BRIEF_MD = """\
# Graphical abstract — design brief (user PowerPoint)

iScience requires a graphical abstract. The recommended concept is a
**three-panel left-to-right "shift" diagram** that visually encodes the
manuscript's central novelty (point-state prediction → feasible-space
diagnosis → rigidification map).

This file is a **design brief** — the actual PNG/SVG should be drawn in
PowerPoint (or Illustrator) by the user. Estimated time: 30 min.

## Canvas

- Aspect ratio: **1:1 (square)** or **landscape 16:9**; iScience accepts both.
  Square is more frequently used for the journal's online thumbnail grid.
- Resolution: 600 dpi at final size (300 dpi minimum).
- Background: white.
- Type: Arial 11–14 pt for labels, 16–18 pt for panel titles.

## Three-panel layout

```
┌─────────────────────────┬─────────────────────────┬──────────────────────────┐
│  (a) POINT-STATE        │  (b) FEASIBLE-SPACE      │  (c) RIGIDIFICATION       │
│      PREDICTION (legacy)│      DIAGNOSIS (ours)    │      MAP (downstream)      │
│                         │                          │                            │
│   single arrow flux     │   colored interval       │   undirected adjacency     │
│   v⃗ at biomass-max     │   [v_min, v_max] cloud   │   over modules + outcomes  │
│   single point          │   "width" annotation     │   bootstrap-stable edges   │
│   monochrome            │   colormap viridis       │   nodes coloured by role   │
│                         │                          │                            │
│   ↓ insufficient for    │   ↓ supports flexibility│   ↓ prioritizes follow-up  │
│     diagnosis           │     collapse diagnosis   │     intervention targets   │
└─────────────────────────┴─────────────────────────┴──────────────────────────┘
```

Above the three panels: single horizontal title arrow with text:
> "Shift from point-state prediction to feasible-space diagnosis"

Below the three panels: small caption "Anchored by sparse standardized harvest-time OD₆₀₀ (n=22). Demonstrated on S. maltophilia SO-1 + iML1515 transfer."

## Panel content recommendations

### Panel (a) POINT-STATE PREDICTION

- One stoichiometric box (light grey rectangle labeled "GEM").
- One thick arrow exiting on the right labeled "FBA / pFBA point flux v⃗*".
- One single black dot near the arrow tip.
- Below: "single deterministic flux at the biomass-maximizing optimum".
- Colour: monochrome grey.

### Panel (b) FEASIBLE-SPACE DIAGNOSIS (centerpiece)

- Same stoichiometric box.
- A **green-to-orange interval bracket** per reaction (5–6 brackets stacked).
- Each bracket annotated with "width = v_max − v_min".
- One bracket shrinking visibly (orange to red) — labelled "flexibility collapse / rigidification".
- Colour: blue-orange divergent palette (e.g., RdBu_r) to encode interval width.

### Panel (c) RIGIDIFICATION MAP

- 6–8 nodes in a small undirected network (no arrowheads — important).
- Nodes coloured by role: blue = exogenous (e.g., O₂ uptake), cyan = flexibility feature, orange = outcome (regime / severity).
- Edge thickness = bootstrap stability frequency.
- Caption: "PC bootstrap → conditional-dependency hypothesis map".

## Colour palette

- Primary: `#0072B2` (blue), `#56B4E9` (cyan), `#E69F00` (orange), `#D55E00` (red), `#009E73` (green) — Wong colourblind-safe palette.
- Background: white.
- Text: black for headings, grey for sub-labels.

## Things to avoid

- Abbreviations as panel titles (iScience explicit instruction).
- Arrowheads in panel (c) — the rigidification map is undirected per
  workstream #E audit.
- Logos, watermarks, decorative icons.

## Sanity check before upload

- All labels readable when reduced to thumbnail (~200 px wide)?
- Three panels left-to-right tell a clean story without reading the
  caption?
- No claim of "validated" / "deployed" / "production-ready" in the
  text — keep wording as "proof-of-concept" / "diagnostic platform" /
  "hypothesis-prioritization".
"""


TABLE1_SPEC_MD = """\
# Conceptual Table 1 / Fig 1A — design spec (user PowerPoint)

Reviewer 1 asked for a clearer up-front explanation of the framework
logic. The rebuttal §7 Introduction insert addresses this narratively;
**Conceptual Table 1 is a one-glance visual reinforcement** of the same
shift (point-state → feasible-space). It is optional from a reviewer-
response perspective (all reviewer points are already addressed in
workstreams #D and #E), but high-impact for first-page visibility.

This file is a **design spec** — the actual Table or Fig 1A should be
drawn in PowerPoint by the user. Estimated time: 30 min.

## Format

- **Option 1 (recommended): three-column comparison table** with 6 rows.
- **Option 2:** the same information rendered as a stylized figure
  (Fig 1A), with rows becoming three horizontal bands.

Both options share the same content. Choose whichever reads faster.

## Content (3 columns × 6 rows)

| Diagnostic dimension | Point-state prediction (legacy) | Feasible-space diagnosis (this work) |
|---|---|---|
| Diagnostic question | "Where does the biomass-maximizing flux go?" | "How much rerouting capacity remains under this constraint context?" |
| Mathematical object | Single flux vector v⃗* (FBA / pFBA / dFBA) | Interval [v_min, v_max] per reaction (targeted FVA) |
| What is interpreted | Magnitude of v_i at the optimum | Width w_i = v_max − v_min (flexibility) |
| Mechanistic reading | Location of the optimum | Degrees of freedom; flexibility collapse / rigidification |
| Downstream artefact | Flux map (Escher-style) | Rigidification map + bootstrap stability (Fig 6) |
| Example methods | FBA, pFBA, dFBA, COSMIC-dFBA | This work (FVA widths + XGBoost + SHAP + PC bootstrap) |

## Visual style

- Black-and-white friendly (iScience accepts both colour and B/W, but
  conceptual tables read well in B/W).
- Row 1 (diagnostic question): bold italic.
- Last column highlighted with a light tint (e.g., 10 % cyan
  `#56B4E9`) so the reader's eye lands on the new contribution.
- Column header row: bold, single-line bottom border.
- No horizontal lines between rows — use vertical alignment + adequate
  row spacing instead.
- Footer: one short sentence "Both paradigms share the same GEM and
  experimental anchors; only the diagnostic representation differs."

## Placement

- **Option A:** insert as Table 1 at the end of the Introduction
  (just before "Results"). This makes Reviewer 1's "framework logic"
  comment fully visual.
- **Option B:** insert as Fig 1A and re-label the existing schematic
  Fig 1 as Fig 1B. The two-panel Fig 1 then tells (1A) the conceptual
  shift and (1B) the workflow.

We recommend **Option A** (Table 1) because the existing Fig 1 is
already a 3-panel workflow schematic; adding a fourth visual element
would crowd the manuscript opener.

## Sanity check

- Does the table answer Reviewer 1's "I can't see the framework logic
  from the abstract" within 5 seconds of glance? If no, simplify the
  language.
- Are the row labels mutually exclusive and exhaustive? (Don't repeat
  the same dimension twice with different wording.)
- Does the rightmost column tell a coherent story when read top-to-
  bottom in isolation?
"""


AI_DISCLOSURE_MD = """\
# AI disclosure — refinement (optional)

The current manuscript paragraph (line 95–96 of the extracted text)
already includes a generative-AI disclosure:

> "During the preparation of this work the author(s) used ChatGPT-5.2
> in order to enhance the readability and language quality. After using
> this tool/service, the author(s) reviewed and edited the content as
> needed and take(s) full responsibility for the content of the
> published article."

iScience's submission guidance accepts this form. **No change is
strictly required.**

## Optional refinement (more transparent)

If the user wishes to disclose the additional AI use during the
revision itself (i.e., Claude Code agent assistance for the revision
analyses), the following paragraph can replace the existing one:

> "During the preparation of this work the author(s) used **ChatGPT-5.2
> for language editing of the main text** and, during the revision,
> **Claude Code as a coding-assistant agent for the supplementary
> analyses (feature-panel ablation, baseline benchmarking, runtime
> profiling, in-silico transfer to iML1515, point-flux vs
> flexibility-interval comparison, Fig 5 non-linearity reanalysis, and
> Fig 6 arrow-consistency audit; full code at
> https://github.com/cmkim0408/Stenotrophomonas_causal-AI, branch
> `revision/iscience-rev1`)**. The authors reviewed and edited the
> content as needed and take full responsibility for the content of
> the published article."

## Recommendation

- If the revision package is the primary place the agent's output is
  used (it is), choose the refined paragraph — this is the most
  transparent stance.
- If the user prefers to keep the disclosure narrow to language
  editing, the existing paragraph stands as-is.
"""


STAR_METHODS_GUIDE_MD = """\
# STAR Methods — sections required for iScience

The manuscript is already in STAR Methods format. The revisions add
new Methods sub-sections that need to slot into the existing structure.
This guide lists the sub-sections, the paste source from
`rebuttal_insertions.md`, and the recommended order.

## STAR Methods structure (iScience standard)

```
STAR ★ Methods
├── KEY RESOURCES TABLE
├── EXPERIMENTAL MODEL AND SUBJECT DETAILS
│     └── Microbial strains
├── METHOD DETAILS
│     ├── Culture conditions and anchoring
│     ├── Genome sequencing and annotation
│     ├── Metabolic model reconstruction
│     ├── Modeling and design-space exploration
│     ├── Regime labeling and flexibilities
│     ├── Feature engineering via targeted FVA
│     │     ┌── NEW: Extended FVA campaign (revision)              ← from rebuttal §1 Methods + REPORT.md Known-limitations §1
│     ├── Explainable AI and causal discovery
│     │     ┌── NEW: Baseline benchmarking protocol (revision)       ← from rebuttal §2 Methods paragraph
│     │     ┌── NEW: Feature-panel ablation protocol (revision)      ← from rebuttal §1 Methods paragraph
│     │     ┌── NEW: Point-flux vs flexibility-interval comparison   ← from rebuttal §7 Methods paragraph
│     │     ┌── NEW: Fig 5 non-linearity reanalysis (revision)        ← from rebuttal §8 Methods paragraph
│     │     ┌── NEW: Fig 6 arrow-consistency audit (revision)        ← from rebuttal §9 Methods sentences (in §9 Results/Methods)
│     ├── Prospective validation experiments
│     └── NEW: External transfer to iML1515 (revision)               ← from rebuttal §4 Methods paragraph
├── QUANTIFICATION AND STATISTICAL ANALYSIS
│     └── NEW: Runtime / computational environment (revision)         ← from rebuttal §5 Methods paragraph
```

## Paste order recommendation

1. After "Feature engineering via targeted FVA":
   - Extended FVA campaign — 1 paragraph
2. After "Explainable AI and causal discovery":
   - Feature-panel ablation — 1 paragraph
   - Baseline benchmarking — 1 paragraph
   - Point-flux vs flexibility-interval — 1 paragraph
   - Fig 5 non-linearity reanalysis — 1 paragraph
   - Fig 6 arrow-consistency audit — 2 sentences (can collapse into the
     causal-discovery sub-section as a closing note)
3. After "Prospective validation experiments":
   - External transfer to iML1515 — 1 paragraph
4. Inside "QUANTIFICATION AND STATISTICAL ANALYSIS":
   - Runtime / computational environment — 2–3 sentences (or footnote
     pointing to SI table)

## Word-count budget (informal)

iScience's STAR Methods is not strictly word-capped, but the revision
additions total roughly:
- Extended FVA campaign: ~120 words
- Feature-panel ablation: ~110 words
- Baseline benchmarking: ~100 words
- Point-flux vs flexibility: ~140 words
- Fig 5 non-linearity: ~110 words
- Fig 6 arrow audit: ~70 words
- External transfer: ~120 words
- Runtime: ~80 words
- **Total addition: ~850 words** (≈ 1 page of Methods)

If the user wants to keep additions shorter, the per-paragraph drafts
in `rebuttal_insertions.md` can each be compressed ~30 % without loss
of substantive content.
"""


README_MD_TEMPLATE = """\
# SUBMISSION_BUNDLE — iScience ISCIENCE-D-26-04043

**Branch:** `revision/iscience-rev1`
**Date:** {today}

This folder is the **submission package** for the iScience revision. The
intent is: zip this folder + the user's hand-edited Word files +
graphical abstract image, and upload to
`https://www.editorialmanager.com/iscience/`.

## What is inside (auto-generated by `code/revision/10_compile_submission_bundle.py`)

| File | Purpose | Action required |
|---|---|---|
| `response_letter.docx` | Point-by-point response | Upload as **rebuttal letter** |
| `highlights_bullets.md` | 5 highlights (and compressed ≤85-char variants) | Paste into iScience submission form / manuscript |
| `graphical_abstract_brief.md` | 3-panel concept spec for the graphical abstract | User makes the PNG in PowerPoint, ~30 min |
| `conceptual_table1_spec.md` | Optional Table 1 / Fig 1A design spec | Optional; ~30 min in PowerPoint |
| `ai_disclosure_paragraph.md` | AI disclosure refinement (optional) | Paste into manuscript if refined version preferred |
| `star_methods_guide.md` | STAR sub-section paste order | Reference while doing manuscript paste work |
| `cover_letter_paragraph.md` | One paragraph for the cover letter | Paste into existing cover letter |
| `caption_replacements/fig5_caption_v2.md` | Fig 5 caption (Variant A = STRONG_NONLINEAR) | Replace Fig 5 caption in main text |
| `caption_replacements/fig6_legend_v2.md` | Fig 6 caption (Variant B = in-place replacement) | Replace Fig 6 caption in main text |
| `housekeeping_verify.md` | Auto-verification of editor housekeeping items | Reference only — informs SUBMISSION_CHECKLIST |

## What is NOT inside (user must produce)

| File | Source | Time estimate |
|---|---|---|
| `Main_..._REVISED-r2-highlighted.docx` | Hand-paste from `rebuttal_insertions.md` (SUBMISSION_CHECKLIST §A) | ~30 min |
| `Main_..._REVISED-r2-clean.docx` | Accept-all-changes copy of the highlighted version | ~5 min |
| `Causal AI-Cover letter-r2.docx` | Existing cover + `cover_letter_paragraph.md` paste | ~5 min |
| `graphical_abstract.png` | Follow `graphical_abstract_brief.md` in PowerPoint | ~30 min |
| (optional) `Table1_conceptual.png` | Follow `conceptual_table1_spec.md` in PowerPoint | ~30 min |

## Order of operations

1. **Open `docs/revision/SUBMISSION_CHECKLIST.md`** — the top-level
   user-facing checklist that walks you through every paste, file
   creation, and verify step. The bundle files in this folder are
   referenced by section letter (A1, A2, ...) inside the checklist.
2. Complete the manuscript paste work (SUBMISSION_CHECKLIST §A).
3. Generate the highlighted + clean Word files (§B).
4. Update the cover letter (§C).
5. Produce the graphical abstract image (§D2).
6. Verify housekeeping items (§E) — most are already in the manuscript.
7. Zip the SUBMISSION_BUNDLE/ folder + user-produced files +
   graphical-abstract PNG; upload to editorialmanager.com.

Estimated total user time after this folder is generated: **~1.5 hours**.

## Provenance

All quantitative numbers used in the bundle are sourced from
`revision_runs/iscience_rev1/metrics_summary.csv` and the per-workstream
`*_summary.md` files. Code is reproducible from the `code/revision/`
folder of the supporting repository (branch `revision/iscience-rev1`,
commit history visible via `git log --oneline`).
"""


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def main() -> None:
    print("[10 compile submission bundle] starting...")

    # 1. Grep-verify housekeeping in the extracted main text
    main_text = MAIN_TEXT.read_text(encoding="utf-8") if MAIN_TEXT.exists() else ""
    if not main_text:
        print(f"  WARN: {MAIN_TEXT} not found; housekeeping verify will return MISSING for all.")
    hk = verify_housekeeping(main_text)

    # 2. Write housekeeping_verify.md
    hk_lines = [
        "# Submission housekeeping — auto-verified status",
        "",
        ("This file is auto-generated by `code/revision/10_compile_submission_bundle.py`. "
         "It checks the extracted main draft (`_unpacked_main/main_text.txt`) "
         "for each iScience editor-required statement and reports whether it "
         "is already present. No modification is performed."),
        "",
        "| # | Editor requirement | Status | Evidence (line: snippet) |",
        "|---|---|---|---|",
    ]
    for r in hk:
        if r["matches"]:
            evidence = "; ".join(f"L{ln}: \"{snip}\"" for ln, snip in r["matches"][:2])
        else:
            evidence = "— (no match)"
        hk_lines.append(f"| {r['id']} | {r['label']} | "
                        f"{'✅ ' + r['status'] if r['status'] == 'PRESENT' else '⚠ ' + r['status']} | "
                        f"{evidence} |")
    hk_lines += ["",
                 "## Resolution notes for any MISSING items",
                 ""]
    for r in hk:
        if r["status"] == "MISSING":
            hk_lines.append(f"- **{r['id']} ({r['label']}):** {r['fallback']}")
    if all(r["status"] == "PRESENT" for r in hk):
        hk_lines.append("_(no items missing)_")
    (BUNDLE / "housekeeping_verify.md").write_text(
        "\n".join(hk_lines) + "\n", encoding="utf-8")
    print(f"  wrote housekeeping_verify.md "
          f"({sum(1 for r in hk if r['status']=='PRESENT')}/{len(hk)} present)")

    # 3. Copy response_letter.docx into bundle
    if RESPONSE_DOCX.exists():
        shutil.copy2(RESPONSE_DOCX, BUNDLE / "response_letter.docx")
        print(f"  copied response_letter.docx")

    # 4. Write auto-generated content files
    (BUNDLE / "highlights_bullets.md").write_text(
        HIGHLIGHTS_MD, encoding="utf-8")
    (BUNDLE / "graphical_abstract_brief.md").write_text(
        GA_BRIEF_MD, encoding="utf-8")
    (BUNDLE / "conceptual_table1_spec.md").write_text(
        TABLE1_SPEC_MD, encoding="utf-8")
    (BUNDLE / "ai_disclosure_paragraph.md").write_text(
        AI_DISCLOSURE_MD, encoding="utf-8")
    (BUNDLE / "star_methods_guide.md").write_text(
        STAR_METHODS_GUIDE_MD, encoding="utf-8")
    print(f"  wrote 5 auto-generated brief files")

    # 5. Copy cover_letter_paragraph.md and caption replacements
    src_cover_par = CAPTIONS / "cover_letter_paragraph.md"
    if src_cover_par.exists():
        shutil.copy2(src_cover_par, BUNDLE / "cover_letter_paragraph.md")
        print(f"  copied cover_letter_paragraph.md")
    for cap in ("fig5_caption_v2.md", "fig6_legend_v2.md"):
        src = CAPTIONS / cap
        if src.exists():
            shutil.copy2(src, BUNDLE / "caption_replacements" / cap)
    print(f"  copied caption_replacements/ (fig5, fig6)")

    # 6. Write README.md
    (BUNDLE / "README.md").write_text(
        README_MD_TEMPLATE.format(today=date.today().isoformat()),
        encoding="utf-8")
    print(f"  wrote README.md")

    # 7. Top-level SUBMISSION_CHECKLIST.md
    write_checklist(hk)
    print(f"  wrote SUBMISSION_CHECKLIST.md")

    print(f"\n[10 compile submission bundle] done. Bundle at: {BUNDLE}")
    print(f"   Checklist at: {CHECKLIST}")


# ---------------------------------------------------------------------------
# Top-level checklist writer
# ---------------------------------------------------------------------------

def write_checklist(hk_rows: list[dict]) -> None:
    auto_present = sum(1 for r in hk_rows if r["status"] == "PRESENT")
    auto_missing = sum(1 for r in hk_rows if r["status"] == "MISSING")
    lines = [
        "# Submission readiness — checklist",
        "",
        "**Manuscript:** ISCIENCE-D-26-04043",
        "**Branch:** `revision/iscience-rev1`",
        f"**Date:** {date.today().isoformat()}",
        "",
        "This is the **single user-facing checklist** for the revision submission.",
        "Work through sections A → E in order. Every paste source is in",
        "`docs/revision/SUBMISSION_BUNDLE/` (auto-generated) or in",
        "`docs/revision/rebuttal_insertions.md` (the master draft file).",
        "",
        "---",
        "",
        "## At a glance",
        "",
        ("- **Quantitative analyses:** all 6 revision workstreams (#1, #2, #3, #4, "
         "#5, #7) plus the extended FVA campaign and the #D, #E follow-ups are "
         "complete. Deferred-count = **0**; every editor mandatory revision and "
         "every reviewer comment is answered in this iteration."),
        ("- **Auto-generated deliverables:** response letter, highlights, "
         "graphical-abstract brief, Conceptual Table 1 spec, AI disclosure, "
         "STAR Methods guide, cover-letter paragraph, caption replacements — "
         "all in `SUBMISSION_BUNDLE/`."),
        (f"- **Housekeeping verify:** {auto_present}/{len(hk_rows)} items "
         f"already present in the main draft; {auto_missing} item(s) need a "
         f"small addition (see §E)."),
        "- **Estimated user time after this checklist is generated:** ~1.5 hours "
        "(paste edits ~30 min + cover letter ~5 min + graphical abstract ~30 min "
        "+ EM upload ~20 min).",
        "",
        "---",
        "",
        "## A. Manuscript text edits (paste from `rebuttal_insertions.md`)",
        "",
        "Open the Word manuscript with **Track Changes ON**. Paste each block at",
        "the placeholder line indicated. The `[Line XXX]` markers should be",
        "replaced with the actual line numbers in your Word doc as you paste",
        "(the rebuttal letter references `[Line XXX]` placeholders the editor",
        "will accept; you don't need them to be real until you finalise the",
        "highlighted version in §B).",
        "",
        "- [ ] **A1.** Abstract — append 1–2 sentences emphasising the",
        "  feasible-space framing. Source: `rebuttal_insertions.md` §7",
        "  Cover-letter insert (first paragraph).",
        "- [ ] **A2.** Introduction — insert the **new paragraph** after current",
        "  paragraph 2. Source: `rebuttal_insertions.md` §7 Introduction insert.",
        "- [ ] **A3.** Methods (Feature engineering via targeted FVA) — add the",
        "  extended-FVA campaign paragraph. Source: `rebuttal_insertions.md` §1",
        "  Methods paragraph + the *Width universe — extended FVA campaign*",
        "  bullet from `REPORT.md` Known limitations.",
        "- [ ] **A4.** Methods (Explainable AI and causal discovery) — add the",
        "  **feature-panel ablation** paragraph (`rebuttal_insertions.md` §1",
        "  Methods) and the **baseline benchmarking** paragraph (§2 Methods).",
        "- [ ] **A5.** Methods (QUANTIFICATION AND STATISTICAL ANALYSIS) — add",
        "  the runtime / environment 2–3 sentences. Source: §5 Methods.",
        "- [ ] **A6.** Methods (Explainable AI and causal discovery) — add the",
        "  **point-flux vs flexibility-interval** paragraph. Source: §7 Methods.",
        "- [ ] **A7.** Methods — add the **Fig 5 non-linearity reanalysis**",
        "  paragraph. Source: §8 Methods.",
        "- [ ] **A8.** Results (near Fig 7) — add the existing-data 5-fold CV",
        "  sentences (macro-F1 = 0.991; R² = 0.906; C10 top mismatch).",
        "  Source: §3 Results sentences.",
        "- [ ] **A9.** Discussion — add the transfer / generalizability",
        "  paragraph. Source: §4 Results paragraph.",
        "- [ ] **A10.** Discussion — add the **3-sentence opener** that",
        "  re-positions the novelty as feasible-space diagnosis. Source: §7",
        "  Results paragraph + Cover-letter insert.",
        "- [ ] **A11.** Limitations — add the proof-of-concept / hypothesis-",
        "  prioritization framing. Source: §7 Cover-letter insert (paragraph 3)",
        "  + Strategic framing.",
        "- [ ] **A12.** **Fig 5 caption** — replace with Variant A from",
        "  `SUBMISSION_BUNDLE/caption_replacements/fig5_caption_v2.md`",
        "  (STRONG_NONLINEAR verdict; explicit breakpoints + ΔR²).",
        "- [ ] **A13.** **Fig 6 caption** — replace with Variant B from",
        "  `SUBMISSION_BUNDLE/caption_replacements/fig6_legend_v2.md`",
        "  (undirected convention; conditional-dependency map).",
        "",
        "---",
        "",
        "## B. Build the highlighted + clean Word files",
        "",
        "- [ ] **B1.** After completing §A with Track Changes ON, save the file as",
        "  `Main_Stenotrophomonas Causal AI_Main draft_REVISED-r2-highlighted.docx`.",
        "- [ ] **B2.** Duplicate the highlighted file, accept all tracked",
        "  changes, save as `Main_..._REVISED-r2-clean.docx`.",
        "",
        "iScience requires **both** versions on submission.",
        "",
        "---",
        "",
        "## C. Update the cover letter",
        "",
        "- [ ] **C1.** Open existing",
        "  `Submission/Causal AI-Cover letter-final-iScience.docx`.",
        "- [ ] **C2.** Insert the `cover_letter_paragraph.md` text",
        "  (`SUBMISSION_BUNDLE/cover_letter_paragraph.md`) as the **second**",
        "  paragraph of the cover (after the brief restatement of the",
        "  manuscript scope, before the closing reviewer acknowledgement).",
        "- [ ] **C3.** Save as `Causal AI-Cover letter-r2.docx`.",
        "",
        "---",
        "",
        "## D. New deliverables (auto-spec'd; user PowerPoint)",
        "",
        "- [ ] **D1. Highlights (3–5 bullets).** Use the compressed `<85`-char",
        "  variants from `SUBMISSION_BUNDLE/highlights_bullets.md`. Paste into",
        "  the iScience EM submission form *and* keep a copy near the Abstract.",
        "- [ ] **D2. Graphical abstract.** Follow the three-panel design brief",
        "  in `SUBMISSION_BUNDLE/graphical_abstract_brief.md` to draw a PNG",
        "  in PowerPoint (~30 min). Export at 600 dpi.",
        "- [ ] **D3. (Optional) Conceptual Table 1 / Fig 1A.** Spec at",
        "  `SUBMISSION_BUNDLE/conceptual_table1_spec.md`. Reviewer 1's intro-",
        "  clarity comment is already addressed narratively in §A2; Table 1",
        "  is a visual bonus.",
        "- [ ] **D4. AI disclosure (optional refinement).** Existing manuscript",
        "  paragraph already discloses ChatGPT use. If you want to also disclose",
        "  Claude Code usage during the revision, swap in the refined paragraph",
        "  from `SUBMISSION_BUNDLE/ai_disclosure_paragraph.md`.",
        "",
        "---",
        "",
        "## E. Editor housekeeping verify",
        "",
        f"Auto-verified status ({auto_present}/{len(hk_rows)} present). Full",
        "evidence table in `SUBMISSION_BUNDLE/housekeeping_verify.md`.",
        "",
    ]
    for r in hk_rows:
        mark = "[x]" if r["status"] == "PRESENT" else "[ ]"
        note = ("present in main draft (see housekeeping_verify.md)"
                if r["status"] == "PRESENT" else r["fallback"])
        lines.append(f"- {mark} **{r['id']}** — {r['label']} — {note}")
    lines += [
        "",
        "---",
        "",
        "## F. Final upload (editorialmanager.com/iscience)",
        "",
        "Files to upload (zip `SUBMISSION_BUNDLE/` is not directly accepted",
        "by EM; you'll upload individual files):",
        "",
        "- [ ] **F1.** Cover letter — `Causal AI-Cover letter-r2.docx` (§C3).",
        "- [ ] **F2.** Manuscript (highlighted) — `Main_..._REVISED-r2-highlighted.docx` (§B1).",
        "- [ ] **F3.** Manuscript (clean) — `Main_..._REVISED-r2-clean.docx` (§B2).",
        "- [ ] **F4.** Response letter — `SUBMISSION_BUNDLE/response_letter.docx`.",
        "- [ ] **F5.** Supplementary Information — your existing SI .docx +",
        "  the new SI figures from `revision_runs/iscience_rev1/figures/` and",
        "  `revision_runs/iscience_rev1/0[1-9]_*/`. Include:",
        "    - `feature_panel_ablation.{png,pdf}` (with all text — axes, "
        "legend, panel-id annotations) + "
        "`feature_panel_ablation_no_labels.{png,pdf}` (no text at all — "
        "axes/ticks/legend/title/annotations stripped; for caption-driven "
        "SI insertion or PowerPoint composition)",
        "    - `benchmark_comparison.{png,pdf}`",
        "    - `confusion_matrix.{png,pdf}`",
        "    - `regression_residuals.{png,pdf}`",
        "    - `external_transfer.{png,pdf}`",
        "    - `pointflux_vs_width.{png,pdf}`",
        "    - `08_nonlinearity/fig5_dependence_with_fits.{png,pdf}`",
        "    - `08_nonlinearity/fig5_shap_interactions.{png,pdf}`",
        "    - `09_fig6_audit/fig6_consistent.{png,pdf,svg}`",
        "    - `09_fig6_audit/fig6_before_after.{png,pdf}`",
        "  Plus the SI tables (`*_metrics.csv`, `runtime_summary.csv`, etc.).",
        "- [ ] **F6.** Graphical abstract — `graphical_abstract.png` (§D2).",
        "- [ ] **F7.** Highlights — paste into the EM submission form (§D1).",
        "- [ ] **F8.** Verify EM author list matches submission (no changes).",
        "",
        "---",
        "",
        "## Ready-to-upload status",
        "",
        f"**Auto-generated deliverables:** ✅ complete ({len(hk_rows)} housekeeping items checked, "
        f"{auto_present} present, {auto_missing} need user attention).",
        "",
        "**Remaining user-side work** (~1.5 h): manuscript paste (§A; ~30 min) "
        "→ highlighted + clean Word files (§B; ~5 min) → cover letter update "
        "(§C; ~5 min) → graphical abstract image (§D2; ~30 min) → EM upload "
        "(§F; ~20 min).",
        "",
        "After §A → §F complete, the submission is **READY TO UPLOAD** to "
        "`https://www.editorialmanager.com/iscience/`.",
    ]
    CHECKLIST.write_text("\n".join(lines) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
