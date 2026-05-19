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
