# Presentation review guide

This file records the checks required before presentation outputs are used in
the manuscript.

## Scientific boundary

- Rho is support-specific magnitude stability.
- The `0.5` reference is provisional.
- Exact support aliases are collapsed by membership in scientific displays.
- Profile coherence, pooled coherence, and selected signed association are
  separate quantities.
- Signed and BH outputs are conditional descriptive summaries.
- Absence and presence below a reference value are distinct states.
- Variance composition is fitted under K=3 and zero within-support covariance.

## Figure checks

### Figure 1

- Preserve every item-level rho value and every absent coordinate-dataset cell.
- Define any across-dataset median over present datasets only.
- Keep catalytic positions as reference markers.

### Figure 2

- Use unique supports only.
- Retain all 924 item-to-first-aggregate profiles.
- Keep zero-mean and fixed-mean V9A regimes separate.
- Label first-crossing values as provisional sensitivity references.

### Figure 3

- Plot configured-group magnitude rho against **profile coherence**.
- Do not relabel profile coherence as pooled coherence.
- Show exact support sizes and sufficient numeric precision.

### Figure 4

- Use only the Phase 3E ranking tables.
- Preserve absolute-effect ranking and the original sign separately.
- Show all three run pairs for N = 5, 10, and 20.

### Figure 5

- Retain all 248 coordinate keys.
- Keep `ABSENT`, `PRESENT_BELOW`, and `CROSSING` visually distinct.
- Keep the strict all-four definition separate from crossing everywhere present.

### Figure 6

- Preserve the 32 configured-group rows and 28/2/2 classification.
- Display the K=3 and zero-covariance qualifications.
- Avoid causal or prescriptive language.

### Supplementary Figure S1

- Show per-dataset item profiles without smoothing or interpolation.
- Label `0.5` as a provisional reference if shown.

### Supplementary Figure S2

- Plot one mark per distinct fitted effect in Panel A.
- Make copied item-row representations explicit.
- Keep BH status visually secondary.
- State that post-selection FDR control has not been established.

## Release checklist

1. `python -m compileall -q .`
2. Fatal-error ruff check.
3. Self-contained tests on supported Python versions.
4. Full integration suite with explicit roots.
5. Two deterministic builds compared byte-for-byte where the format permits.
6. PNG, PDF, and SVG inspected at manuscript scale.
7. Caption drafts and numeric audit files regenerated with the figures.
8. Manuscript image, caption, and source commit verified together.
