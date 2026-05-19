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
