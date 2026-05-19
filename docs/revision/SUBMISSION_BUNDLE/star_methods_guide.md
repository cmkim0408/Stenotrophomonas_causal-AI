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
