# STRIDE Figure Suite Review

Review target: `paper_figures` only. Frozen STRIDE outputs in
`/Users/medhasharma/Downloads/forcodexstride/strideouts` were read but not modified. This
pass changed presentation, layout, labels, emphasis, and documentation only.

## Narrative Audit

The main suite now follows a defensible sequence:

1. Figure 1 asks where residue-level reproducibility lies along NS2B-NS3.
2. Figure 2 asks at what spatial scale catalytic-region reproducibility saturates.
3. Figure 3 asks how catalytic domains differ from the rest of the protein.
4. Figure 4 asks where statistically robust signed energetic effects occur.
5. Figure 5 asks how reproducible features are conserved across serotypes.
6. Figure 6 asks what limits reproducibility.

The story is intentionally more restrained than an idealized catalytic-core narrative.
The data support heterogeneity, catalytic coherence, broad signed effects, cross-serotype
conservation structure, and tau2-dominated limits. They do not support claims that
catalytic regions are the highest-rho regions or that significant effects concentrate
there.

## Figure-by-Figure Review

### Figure 1 — Residue-level reproducibility is distributed across NS2B-NS3

Biological question: Where are reproducible dynamics located?

Strengths: Opens with the real sequence-wide rho landscape; chain bands orient the reader;
catalytic triad residues are visible but no longer overclaimed.

Weaknesses: A pure landscape remains partly descriptive. Domain labels are necessarily
compressed because several named NS3 domains are short or adjacent.

Changes made: Shortened domain-track labels, staggered leader lines, retained neutral
domain blocks, and preserved catalytic triad markers as landmarks rather than as a
continuous catalytic-region block.

Remaining weaknesses: The figure answers "where" honestly, but the answer is distributed
rather than narratively punchy. The broad "unassigned" domain span remains biologically
unsatisfying.

Publication readiness: 7.5/10.

### Figure 2 — Catalytic-region rho saturates above the residue scale

Biological question: At what spatial scale does reproducibility emerge?

Strengths: The plateau is now visually obvious. A shaded SS/motif/domain plateau region
makes clear that these labels collapse to the same aggregated signal in this hierarchy.
Panel B correctly shows that gates remain mostly residue or secondary-structure.

Weaknesses: The figure needs caption discipline because SS/motif/domain are structurally
present but not independently informative here.

Changes made: Added a plateau band and explanatory inline annotation; softened individual
trajectories; retitled Panel A around the residue-to-aggregate jump.

Remaining weaknesses: The figure is strong, but the scale labels could still be misread
without a caption noting empty SS/motif biology.

Publication readiness: 8.5/10.

### Figure 3 — Catalytic domains are coherent but not the highest-rho regions

Biological question: How do catalytic regions differ from the remainder of the protein?

Strengths: Much more coherent conceptually. Panel A establishes that catalytic domains are
not top-ranked by rho; Panel B now uses grey background points and orange catalytic points
so the reader sees the actual distinction: coherence, not rho.

Weaknesses: This remains conceptually subtle because the result is mixed. "Unassigned"
and Gly45 regions dominate rho, while catalytic domains stand out mainly in coherence.

Changes made: Reframed Panel B as catalytic-vs-other rather than a second serotype-coded
scatter; reduced competing colours; labelled catalytic domain groups directly; shortened
the Panel B title.

Remaining weaknesses: Because the domain heatmap and coherence scatter use different
geometries, they still require the reader to integrate two encodings.

Publication readiness: 8/10.

### Figure 4 — Significant signed effects are widespread across NS2B-NS3

Biological question: Where do statistically robust energetic effects occur?

Strengths: The visual hierarchy now matches the title. FDR-significant points dominate in
black; non-significant points recede; catalytic highlighting is contextual rather than
overpowering.

Weaknesses: A volcano plot is efficient but abstract; positional context is left to the
supplementary figure.

Changes made: Restricted orange outlines to catalytic rows that are also FDR-significant;
reduced outline size/alpha; changed labels to the strongest FDR-significant effects only.

Remaining weaknesses: Dense regions on the positive-effect side still contain many points,
but this density is the data, not a plotting artifact.

Publication readiness: 8/10.

### Figure 5 — Cross-serotype reproducibility is heterogeneous across positions

Biological question: How conserved are reproducible features across serotypes?

Strengths: The figure now stays on one question. Panel A gives counts by number of
serotypes; Panel B honestly shows weak profile correlations; Panel C shows that positions
reproducible in more serotypes have higher median residue-level rho.

Weaknesses: Panel B is a qualifying result rather than a headline result. It is useful
because it prevents overclaiming, but it weakens the emotional arc of the figure.

Changes made: Replaced the old signed-effect Panel C with a conservation-class rho
summary. Removed the S4 signed-effect dependency from Figure 5. Kept catalytic triad
markers as landmarks in Panels A and C.

Remaining weaknesses: The figure is now honest and readable, but less dramatic. If space
is tight, Panel B could move to supplementary material.

Publication readiness: 8/10.

### Figure 6 — Reproducibility is limited by between-replicate variance

Biological question: What fundamentally limits reproducibility?

Strengths: The cleanest mechanistic figure. Tau2 dominance is immediately visible across
serotypes and regions.

Weaknesses: Four small multiples repeat y-axis labels and show a 100% stacked bar, which
double-encodes one degree of freedom.

Changes made: Removed star prefixes from catalytic y-axis labels and kept catalytic
emphasis to tick colour only.

Remaining weaknesses: A single mean frac_tau2 panel would use less ink, but the current
facets preserve per-serotype transparency.

Publication readiness: 8.5/10.

### Supplementary Figure S1 — Per-serotype residue-level rho landscapes are heterogeneous

Strengths: Useful backup for Figure 1; shows the raw per-serotype texture that the main
median landscape compresses.

Weaknesses: Still visually busy by design, because it displays all per-serotype residue
profiles.

Changes made: Reduced catalytic-ring weight/alpha and kept the title descriptive.

Remaining weaknesses: Best treated as a supplemental audit, not a main narrative figure.

Publication readiness: 7/10.

### Supplementary Figure S2 — Signed dynamic effects are distributed along NS2B-NS3

Strengths: Now complements Figure 4 without contradicting it. It shows positional spread
while Figure 4 shows significance and magnitude.

Weaknesses: The NS3 signed-effect field is dense; some overlap is unavoidable at this
width.

Changes made: Restricted catalytic outlines to FDR-significant catalytic rows and removed
catalytic status from the label-selection rule.

Remaining weaknesses: Still more crowded than an ideal supplement; a faceted chain/domain
version could be cleaner but would be a larger redesign.

Publication readiness: 7/10.

## Reviewer-Style Comments

The revised suite is substantially more credible. The authors have stopped forcing a
catalytic-localization story onto data that are plainly heterogeneous. This is essential:
in the previous versions, Figures 1, 3, 4 and 5 each risked immediate reviewer rejection
because the title made a stronger claim than the plotted data.

The strongest main-text figures are Figure 2 and Figure 6. They answer clear questions,
their visual encodings are simple, and their conclusions are directly visible. Figure 3 is
also improved because it now communicates a mixed but interesting result: catalytic domains
are not the most reproducible by rho, but they are relatively coherent in signed direction.

The remaining strategic weakness is that the biological narrative is not a single clean
arc toward a conserved catalytic core. Instead, the honest story is: reproducibility is
distributed; scale aggregation matters; catalytic regions differ more by coherence than by
rho; robust signed effects are broad; conservation is heterogeneous but associated with
higher median rho; replicate disagreement limits reproducibility. That is publishable, but
the manuscript must embrace this complexity.

For a Nature/Science-style review, I would still ask for explicit caption caveats:
rho-star is provisional, residue-level claims are exploratory at K=3, SS/motif levels are
not independently biologically annotated in this hierarchy, and serotype-level comparisons
are descriptive with n=4 biological systems.

## Future Work

- Rewrite manuscript Results text to match the revised figure sequence.
- Consider moving Figure 5B to supplementary if the main narrative needs a tighter visual
  arc.
- Consider a future upstream annotation pass to reduce the dominance of "unassigned"
  structural regions. This should be a new analysis/annotation step, not a plotting patch.
- At final typesetting, check all panels at actual journal column widths and confirm label
  size remains >= 6 pt.
- If the manuscript needs a stronger biological hook, develop it around coherence and
  replicate-variance limits rather than unsupported catalytic rho peaks.

## Overall Assessment

Current strongest-to-weakest ranking:

1. Figure 2
2. Figure 6
3. Figure 5
4. Figure 3
5. Figure 4
6. Figure 1
7. Supplementary Figure S1
8. Supplementary Figure S2

Overall publication readiness: 8/10 for a figure suite, assuming captions and manuscript
claims are rewritten to match. The figures are now scientifically honest, readable, and
mostly publication-quality. The main remaining barrier is narrative framing, not plotting
mechanics.
