# Fig 5 — revised caption (workstream #D)

**RECOMMENDED based on workstream #D verdict: `STRONG_NONLINEAR`** → use **Variant A** below.

The non-linearity claim in the published manuscript is *supported* by the
reanalysis. Both the piecewise-linear (knot grid-searched between the 10th
and 90th percentile of each feature) and degree-2 polynomial fits
significantly improve over the linear baseline on the top-3 OOF SHAP
features (p < 1e-4 in every comparison). The dominant signal is in
`width__ADCS` (Δ R² piecewise vs linear = +0.425; breakpoint = 0.00),
followed by `width__ACONT` (Δ R² = +0.053; breakpoint = 1973.79) and
`width__12DGR120tipp` (Δ R² = +0.166; breakpoint = 3.06). The ADCS
breakpoint at zero is a data-driven *flexibility-collapse tipping point*
— exactly the rigidification interpretation the manuscript advances.

Three variants are provided. Variant A is the recommended replacement for
the published Fig 5b/5c caption; Variants B and C are retained for
completeness in case the user prefers a more conservative or
interaction-only framing.

---

## Variant A — STRONG_NONLINEAR / THRESHOLD_LIKE (RECOMMENDED)

> **Figure 5. Determinants of growth-potential variation and data-driven
> non-linear tipping points.** (A) Global importance ranking of FVA-width
> features derived from the XGBoost regressor trained on the normalized
> growth-potential index (severity = obj/obj_max), ranked by mean
> out-of-fold |SHAP| under 5-fold cross-validation (random_state = 42;
> Methods). (B–C) SHAP dependence for the two top-ranked features
> (`width__ADCS`, `width__APSR` in the original 120-width universe; on
> the extended ~300-width universe (workstream #B) the top-ranked
> features become `width__ACONT`, `width__12DGR120tipp`, `width__ADCS`,
> which span the same TCA / acetate-uptake / amino-acid biosynthesis
> modules). Solid orange curve: LOESS overlay (frac = 0.4).
> Dashed green: continuous piecewise-linear fit with the breakpoint
> (vertical dotted line) chosen by grid search; F-test piecewise vs
> single linear: p < 1e-4 in every case. Dot-dashed pink: degree-2
> polynomial fit, F-test polynomial vs linear: p < 1e-3 in every case.
> The piecewise break for `width__ADCS` at width ≈ 0 corresponds to a
> data-driven **flexibility-collapse tipping point** — the SHAP slope
> changes sharply when the feasible interval of ADC synthase contracts
> to zero, consistent with the manuscript's rigidification interpretation
> (Fig 6). See `revision_runs/iscience_rev1/08_nonlinearity/` for the
> full per-feature test table (`nonlinearity_metrics.csv`) and SHAP
> interaction analysis.

### Suggested main-text wording (replaces "non-linear effects" with quantitatively supported language)

> "...SHAP dependence reveals **data-driven tipping points**, with
> piecewise-linear breakpoints (p < 1e-4 vs linear baseline) at
> `width__ADCS` ≈ 0 (flexibility-collapse threshold; Δ R² = +0.425),
> `width__12DGR120tipp` ≈ 3.06 (Δ R² = +0.166), and `width__ACONT` ≈
> 1973.8 (Δ R² = +0.053). The breakpoint at `width__ADCS` ≈ 0 is
> consistent with rigidification setting in as the feasible interval of
> ADC synthase contracts to zero..." (Main text [Line XXX]; precise
> wording at user discretion.)

---

## Variant B — INTERACTION_MOD (not the current verdict; retained for completeness)

> **Figure 5. Determinants of growth-potential variation and
> interaction-modulated dependence.** (A) Global importance ranking of
> FVA-width features derived from the XGBoost regressor (mean OOF |SHAP|,
> 5-fold CV). (B–C) SHAP dependence for the top-ranked features. The
> main dependence appears approximately linear, but stratifies sharply
> when conditioned on the strongest interacting partner reaction (panel
> colour); the average pairwise SHAP-interaction magnitude is reported
> in `revision_runs/iscience_rev1/08_nonlinearity/shap_interactions.csv`.
> This is consistent with metabolic rigidification emerging only when
> both upstream and downstream degrees of freedom contract simultaneously.

---

## Variant C — LINEAR_RETREAT (not the current verdict; retained for completeness)

> **Figure 5. Determinants of growth-potential variation: monotonic
> trends across the LHS condition envelope.** (A) Global importance
> ranking of FVA-width features (mean OOF |SHAP|, 5-fold CV). (B–C) SHAP
> dependence for the top-ranked features. Within the sampled LHS range,
> the dependence is approximately linear (piecewise and polynomial
> F-tests vs linear baseline: p > 0.05); stronger non-linearity may
> emerge at more extreme constraint regimes not sampled by the current
> LHS. We revise the prior "non-linear effects" wording to **"monotonic
> trends consistent with progressive rigidification within the tested
> condition envelope"**.

---

## Note for the response letter

The two-paragraph rebuttal in §8 of `docs/revision/rebuttal_insertions.md`
explicitly cites the verdict (`STRONG_NONLINEAR`), the per-feature
piecewise / polynomial F-test p-values, and the breakpoints. It also
states transparently that the test was performed on the extended
~300-width universe (so the top-ranked features differ slightly from
the published Fig 5b/5c but cover the same TCA / acetate-uptake /
amino-acid-biosynthesis module structure).

The same response closes editor mandatory revision **M8** and reviewer
2's comment **R2.3**.
