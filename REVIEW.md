# STRIDE Figure Suite Review

This review was performed against the real frozen outputs in
`/Users/medhasharma/Downloads/forcodexstride/strideouts`. No STRIDE statistics,
thresholds, parquet files, CSV files, or files under `strideouts` were modified.

## 1. Figure-by-Figure Review

### Figure 1 — Residue-level reproducibility is distributed across NS2B-NS3

Strengths: strong opening orientation; residue-level rho is shown directly; catalytic
triad residues are visible without being overclaimed.

Weaknesses: the previous title claimed catalytic localization, but catalytic-triad
median rho is lower than the non-catalytic distribution. The old domain ribbon also
looked like one continuous catalytic block.

Changes made: retitled the figure; changed the domain track to neutral blocks; kept
catalytic residues as rings/ticks; moved domain labels outside the blocks with leader
lines.

Remaining limitations: the domain track is still compact at double-column width and
uses broad framework domain labels that include "unassigned" spans.

### Figure 2 — Catalytic-region rho saturates above the residue scale

Strengths: the saturation between residue and higher levels is now visually obvious.
Panel B makes the actual gating distribution visible.

Weaknesses: the previous title and subtitle claimed domain-scale emergence, but the
census shows gates remain mostly residue or secondary-structure, with only small chain
counts and no domain-gated class in this table.

Changes made: retitled the figure and both panel subtitles to distinguish trajectory
saturation from gated-resolution counts.

Remaining limitations: secondary-structure and motif levels are structurally present but
biologically empty in this dataset, so the plateau should be described carefully in the
caption.

### Figure 3 — Catalytic domains are coherent but not the highest-rho regions

Strengths: this is now the most scientifically honest main-domain figure. The heatmap
shows that broad unassigned/Gly45 regions have higher rho, while the scatter preserves
the real catalytic coherence signal.

Weaknesses: the old title was unsupported: Catalytic Triad is the lowest mean-rho row
and Oxyanion Loop is mid-ranked.

Changes made: retitled the figure; changed panel titles; preserved catalytic emphasis
only for coherence/functional reference.

Remaining limitations: "unassigned" regions dominate the heatmap and are biologically
less satisfying than named structural modules.

### Figure 4 — Significant signed effects are widespread across NS2B-NS3

Strengths: the volcano remains an efficient summary of signed effect size and FDR
status. Catalytic positions are visible without claiming exclusivity.

Weaknesses: the old title claimed concentration in catalytic regions, but significant
effects are widespread: most FDR-significant signed rows are outside catalytic domains.

Changes made: retitled the figure; reduced the global label cap from 8 to 4 to avoid
label clutter.

Remaining limitations: catalytic rings remain visually dense because many signed rows
fall within catalytic-domain annotations.

### Figure 5 — Cross-serotype reproducibility is heterogeneous across positions

Strengths: now honest and readable. Panel A summarizes all positions without unreadable
y-axis labels; Panel B is kept as a constraint showing weak whole-profile correlations;
Panel C shows signed effects across conservation classes without implying a strong trend.

Weaknesses: the old headline overclaimed a conserved dynamically important core. The
observed conservation-effect association is weak, and serotype profile correlations are
near zero or slightly negative off-diagonal.

Changes made: replaced Panel A's 248 labeled rows with conservation-class counts and
catalytic-triad markers; retitled Panel B; replaced LOWESS with class medians; bounded
deterministic jitter to [0, 1]; reduced label clutter.

Remaining limitations: Panel B deserves inclusion only as negative/qualifying evidence,
not as a headline similarity result. Panel C is descriptive and should not be framed as
a statistical test.

### Figure 6 — Reproducibility is limited by between-replicate variance

Strengths: the core claim is supported: tau2 dominates the variance fraction in most
regions, with median frac_tau2 about 0.93. The panel is clean and interpretable.

Weaknesses: repeated y-axis labels consume space, and the 100% stacked bars double-encode
one degree of freedom because frac_tau2 + frac_sigma2 = 1.

Changes made: no structural code change was required; the title remains supported.

Remaining limitations: a single mean frac_tau2 dot/bar panel would be more compact, but
the current decomposition is clearer for readers.

### Supplementary Figure S1 — Per-serotype residue-level rho landscapes are heterogeneous

Strengths: useful backup for Figure 1; per-serotype traces show real heterogeneity and
avoid cross-serotype aggregation hiding.

Weaknesses: the old description implied catalytic peaks in every serotype, which is not
supported by the traces.

Changes made: retitled and reframed catalytic residues as landmarks.

Remaining limitations: the rolling mean is a visual aid only and should be described as
such if mentioned.

### Supplementary Figure S2 — Signed dynamic effects are distributed along NS2B-NS3

Strengths: useful positional companion to Figure 4; shows that signed effects occur in
both chains and across NS3 regions.

Weaknesses: the old message claimed clustering at catalytic residues, but the significant
effects are broadly distributed.

Changes made: retitled and reduced labels via the shared label rule.

Remaining limitations: dense catalytic rings in the NS3 region can still obscure nearby
points at small print size.

## 2. Scientific Audit of Titles

- Figure 1, "Residue-level reproducibility is distributed across NS2B-NS3": fully
  supported. Top residue-level median-rho positions are spread across NS2B and NS3, and
  catalytic-triad residues are not the maxima.
- Figure 2, "Catalytic-region rho saturates above the residue scale": fully supported.
  F7 trajectories rise from residue and then plateau at SS/motif/domain for this
  hierarchy.
- Figure 3, "Catalytic domains are coherent but not the highest-rho regions": fully
  supported. Catalytic domains have moderate rho but high coherence relative to many
  non-catalytic regions.
- Figure 4, "Significant signed effects are widespread across NS2B-NS3": fully supported.
  FDR-significant signed effects occur across NS2B unassigned, NS3 unassigned, Gly45 Turn,
  Oxyanion Loop, C-Terminal Tail, and Catalytic Triad.
- Figure 5, "Cross-serotype reproducibility is heterogeneous across positions": fully
  supported. Position counts span 0-4 reproducing serotypes and catalytic triad residues
  occupy different classes.
- Figure 6, "Reproducibility is limited by between-replicate variance": fully supported.
  tau2 is the dominant variance fraction in most regions and overall.
- Supplementary Figure S1, "Per-serotype residue-level rho landscapes are heterogeneous":
  fully supported by the per-serotype traces.
- Supplementary Figure S2, "Signed dynamic effects are distributed along NS2B-NS3": fully
  supported by the Manhattan/lollipop distribution.

## 3. Nature Communications Reviewer Critique

The revised suite is materially more honest than the previous render. The strongest
improvement is the removal of catalytic-core overclaiming: several original titles would
have been challenged immediately because the plotted data contradicted them. The current
figures now distinguish reproducibility, coherence, signed significance, and conservation
instead of forcing them into a single catalytic narrative.

The main weakness is that the biological story is less crisp than the original figure
sequence wanted it to be. Figure 5 no longer provides a high-impact "conserved dynamic
core" headline; it provides a sober cross-serotype heterogeneity audit. That is the right
scientific decision, but the manuscript text will need to be rewritten around it. Figure 3
also shows that catalytic domains are not the most reproducible by rho, so any abstract or
Results language saying otherwise should be removed.

For a Nature Communications submission, I would ask for stronger caption discipline:
state that rho* is provisional, that residue-scale claims are exploratory at K=3, and that
SS/motif labels are not biologically populated in this dataset. I would also ask the
authors to justify the inclusion of Figure 5B as a negative result or move it to the
supplement if the manuscript cannot use it interpretively.

## 4. Remaining Work

- Rewrite manuscript captions and Results text to match the revised, weaker but supported
  claims.
- Consider simplifying Figure 6 to a single tau2-fraction summary if journal space is
  tight.
- Check all figures at final journal column sizes after typesetting.
- Decide whether Figure 5B belongs in the main figure as a qualifying panel or should move
  to supplementary material.
- If structural interpretation is central, replace "unassigned" domain labels upstream in
  a future analysis stage; do not patch this in plotting code.

## 5. Overall Assessment

Strongest to weakest:

1. Figure 6
2. Figure 3
3. Figure 2
4. Figure 4
5. Supplementary Figure S1
6. Figure 1
7. Supplementary Figure S2
8. Figure 5

The complete figure suite is closer to publication quality but is not yet fully suitable
for a high-impact submission as a polished narrative package. It is scientifically honest
and technically clean, but the central cross-serotype/catalytic story is weaker than the
previous titles claimed. Submission readiness now depends less on plotting mechanics and
more on rewriting the manuscript argument around heterogeneity, coherence, and variance
limits rather than catalytic localization.
