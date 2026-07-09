# Manuscript Number Audit

Generated from frozen STRIDE outputs under `/Users/medhasharma/Downloads/forcodexstride/strideouts`.

No STRIDE outputs were modified. The long-form source table is `numbers_audit.csv`.

## Source Hierarchy

- Use stage outputs (`outputs_s2`-`outputs_s6`) when a manuscript claim needs columns not exported in S7.
- Use S7 figure/table artifacts when verifying the exact figure/table-ready values.
- For Figure 4 signed effects, the authoritative plotted denominator is `outputs_s4/significance_screen.parquet`, not the provisional `outputs_s2/signed_screen.parquet`.

## Part A - Current Manuscript Claims

### Table 1 / Comparable Substrate

|serotype|n_loci|n_gated_residue|frac_mixed|median_rho|IQR|
|---|---|---|---|---|---|
|DENV1|236|122|0.703|0.515|0.367-0.732|
|DENV2|247|137|0.680|0.535|0.377-0.691|
|DENV3|242|145|0.612|0.558|0.387-0.706|
|DENV4|199|106|0.704|0.528|0.315-0.702|

Status: confirmed. The manuscript values 0.515 / 0.535 / 0.558 / 0.528, n_loci 199-247, n_gated_residue 106-145, and frac_mixed roughly 0.61-0.70 are supported.

### Domain Reproducibility

|domain|rho_range|rho_mean|coherence_range|
|---|---|---|---|
|Catalytic Triad|0.362-0.778|0.526|0.572-0.965|
|Oxyanion Loop|0.565-0.894|0.680|0.712-0.889|
|Unassigned NS3|0.903-0.978|0.958||

Status: confirmed, with one naming clarification. The manuscript phrase "Unassigned NS3" corresponds to `chain == NS3` and `domain == unassigned` in the frozen table.

### Signed Effects / Figure 4

- Per-locus FDR-controlled significance screen: 152/924 loci are FDR-significant.
- Signed loci in that table: 302/924.
- Catalytic Triad FDR-significant loci: 3/152.
- Oxyanion Loop FDR-significant loci: 11/152.
- Catalytic Triad + Oxyanion Loop total: 14/152 (9.2% of FDR-significant loci).
- Provisional S2 screen-pass / CI-excluding-zero count at rho*=0.5: 181. This is not the Figure 4 FDR denominator.

|chain|domain|loci|signed|fdr|
|---|---|---|---|---|
|NS3|unassigned|524|141|78|
|NS2B|unassigned|203|65|34|
|NS3|Gly45 Turn|32|24|22|
|NS3|Oxyanion Loop|32|28|11|
|NS3|C-Terminal Tail|85|35|4|
|NS3|Catalytic Triad|12|5|3|
|NS3|120s Loop|24|3|0|
|NS3|B2b-C2 Hairpin|12|1|0|

Interpretation: "significant signed effects are widespread" is supported. "Catalytic minority" is supported by absolute FDR counts under the Figure 4 denominator (14/152); however, it should not be worded as a formal enrichment/depletion claim. Descriptively, catalytic+oxyanion regions have 14/44 FDR-significant loci (31.8%) versus 138/880 outside those domains (15.7%), so the outside majority is driven by many more outside loci, not by a lower catalytic per-locus rate.

### Cross-Serotype Conservation / Figure 5

- Conservation-class counts for 0, 1, 2, 3, 4 reproducing serotypes: 18, 49, 100, 63, 18.
- Positions reproducible in all four serotypes: 18.
- Catalytic-triad residues among all-four positions: 0.
- Median rho by conservation class 0-4: 0.339, 0.395, 0.529, 0.625, 0.744.
- Monotonic median-rho increase across classes: True.

All-four-serotype positions:

|position|chain|domain|median_rho|signed_serotypes|direction|catalytic|
|---|---|---|---|---|---|---|
|NS2B:-17|NS2B|unassigned|0.754|2|none|False|
|NS2B:-19|NS2B|unassigned|0.709|2|none|False|
|NS2B:-23|NS2B|unassigned|0.720|2|none|False|
|NS2B:-24|NS2B|unassigned|0.770|1|NA|False|
|NS2B:-25|NS2B|unassigned|0.617|2|none|False|
|NS2B:-9|NS2B|unassigned|0.745|2|none|False|
|NS3:121|NS3|120s Loop|0.625|1|NA|False|
|NS3:43|NS3|Gly45 Turn|0.648|2|decrease|False|
|NS3:47|NS3|Gly45 Turn|0.771|3|decrease|False|
|NS3:106|NS3|unassigned|0.767|2|increase|False|
|NS3:130|NS3|unassigned|0.772|2|none|False|
|NS3:131|NS3|unassigned|0.743|1|NA|False|
|NS3:33|NS3|unassigned|0.677|2|none|False|
|NS3:34|NS3|unassigned|0.666|1|NA|False|
|NS3:72|NS3|unassigned|0.788|3|increase|False|
|NS3:84|NS3|unassigned|0.727|2|none|False|
|NS3:85|NS3|unassigned|0.789|2|none|False|
|NS3:94|NS3|unassigned|0.794|2|none|False|

Status: confirmed. The all-four-serotype core is limited and non-catalytic in the frozen table.

### Variance / Figure 6

- Median frac_tau2: 0.925.
- Variance regimes: replicate_dominated=28, sampling_dominated=2, balanced=2.
- Overall frac_tau2 range: 0.358-0.994.
- Catalytic Triad frac_tau2 range: 0.834-0.968.
- Oxyanion Loop frac_tau2 range: 0.703-0.981.

Non-replicate-dominated exceptions:

|serotype|chain|domain|frac_tau2|frac_sigma2|regime|
|---|---|---|---|---|---|
|DENV1|NS3|120s Loop|0.361|0.639|sampling_dominated|
|DENV1|NS3|Gly45 Turn|0.594|0.406|balanced|
|DENV3|NS3|Gly45 Turn|0.425|0.575|balanced|
|DENV4|NS3|B2b-C2 Hairpin|0.358|0.642|sampling_dominated|

Status: median frac_tau2 is confirmed, but the manuscript's "29/32 replicate-dominated" and "three sampling-dominated exceptions" wording is not supported by the stored regime labels. The frozen table says 28 replicate-dominated, 2 balanced and 2 sampling-dominated.

### Hierarchy / Figure 2

Provisional Figure 2B gated-scale categories: chain, residue, secondary_structure. Domain is not present as a provisional census category.

|serotype|scale|n_loci|
|---|---|---|
|DENV1|residue|122|
|DENV1|secondary_structure|114|
|DENV2|residue|137|
|DENV2|secondary_structure|105|
|DENV2|chain|5|
|DENV3|residue|145|
|DENV3|secondary_structure|94|
|DENV3|chain|3|
|DENV4|residue|106|
|DENV4|secondary_structure|86|
|DENV4|chain|7|

Scale-curve rho summaries:

|scale|mean_rho|median_rho|range|n|
|---|---|---|---|---|
|chain|0.947|0.973|0.645-0.990|924|
|complex|0.963|0.993|0.863-0.993|924|
|domain|0.897|0.975|0.344-0.979|924|
|motif|0.897|0.975|0.344-0.979|924|
|protein|0.963|0.993|0.863-0.993|924|
|residue|0.527|0.537|0.014-0.936|924|
|secondary_structure|0.897|0.975|0.344-0.979|924|

Status: the plateau claim is supported. Secondary-structure, motif and domain have identical mean/median/range summaries in the full scale curve, and catalytic-region F7 has the same collapse across those three aggregate levels.

### Provenance

- K=3 per serotype: confirmed from `outputs_s6/replicate_regime.parquet`.
- Provisional rho*: 0.5.
- Serotypes present: DENV1, DENV2, DENV3, DENV4.
- S7 artifacts: 13 total (8 figures and 5 tables).

## Part B - Additional Useful Numbers

### Dataset / Case-Study Scale

- Serotypes: 4.
- Replicate runs: 12 total, K=3 per serotype.
- Replicate-table rows: 2772.
- Serotype-residue rows: 924; union residue positions: 248; positions present in all serotypes: 199.
- Union positions by chain: NS2B=62, NS3=186.
- Domain-by-serotype entries: 32.
- Cross-serotype position rows: 248.
- Per-locus signed effects: 302; mixed mechanisms at provisional rho*: 622; FDR-significant signed effects: 152.

### Reproducibility Landscape

- All residue-level rho values: median 0.537, IQR 0.368-0.710, range 0.014-0.936.
- Chain medians: NS2B=0.559, NS3=0.527.
- Top/bottom 10 median-rho positions and catalytic-triad residue-level rho values are enumerated in `numbers_audit.csv`.

### Spatial Hierarchy

- Mean rho gain from residue to first aggregate: 0.369.
- Residue mean/median rho: 0.527/0.537.
- Secondary-structure/motif/domain mean/median rho: 0.897/0.975.

### Domain Analysis

|rank|chain|domain|mean_rho|rho_range|mean_coherence|
|---|---|---|---|---|---|
|1|NS3|unassigned|0.958|0.903-0.978|0.234|
|2|NS3|Gly45 Turn|0.885|0.768-0.979|0.661|
|3|NS2B|unassigned|0.854|0.645-0.978|0.457|
|4|NS3|120s Loop|0.777|0.517-0.972|0.932|
|5|NS3|C-Terminal Tail|0.713|0.495-0.939|0.489|
|6|NS3|Oxyanion Loop|0.680|0.565-0.894|0.792|
|7|NS3|B2b-C2 Hairpin|0.615|0.344-0.837|0.767|
|8|NS3|Catalytic Triad|0.526|0.362-0.778|0.760|

Useful contrast: NS3 unassigned has the highest mean rho but low coherence; Catalytic Triad has the lowest mean rho among listed domains but relatively high coherence.

### Signed Effects

- FDR-significant effects by sign: negative=73, positive=79.
- Top positive/negative beta loci, top adjusted-p loci, and catalytic-domain signed beta values are enumerated in `numbers_audit.csv`.

### Cross-Serotype Conservation

- Exactly-three-serotype reproducible positions: 63; all are listed in `numbers_audit.csv`.
- Serotype/chain reproducible-position counts are listed in `numbers_audit.csv`.

### Variance Budget

- Overall frac_tau2 IQR: 0.833-0.963.
- Replicate-dominated fraction: 28/32 (87.5%).
- Sampling-dominated fraction: 2/32 (6.2%).
- Most replicate- and sampling-dominated regions are enumerated in `numbers_audit.csv`.

## Missing or Not Verifiable from Frozen Outputs

- Sequence identity 65-79% and static-structure RMSD below 1.5 A are background claims and are not verifiable from the frozen STRIDE reduction outputs.
- MD engine/version, force field, water model, ensemble, temperature/pressure/ions, equilibration, trajectory length, frame sampling, PDB accessions and shared-reference-frame control results are not present in the inspected frozen reduction outputs.
- The manuscript's Figure 5 caption still says catalytic-triad residues are marked in Panel A; after the presentation cleanup this should be edited if the figure no longer marks catalytic residues there.
