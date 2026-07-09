# Numbers For Manuscript

## Results 4.1

- Four serotypes, K=3 independent replicate runs per serotype (12 total replicate runs).
- Analysable loci per serotype: 236 / 247 / 242 / 199 for DENV1-DENV4 (range 199-247).
- Gated residues at provisional ρ*=0.5: 122 / 137 / 145 / 106 (range 106-145).
- Median residue-level ρ: 0.515 / 0.535 / 0.558 / 0.528 for DENV1-DENV4.
- Mixed-mechanism fraction: 0.703 / 0.680 / 0.612 / 0.704.

## Figure 1 / Landscape

- All residue-level ρ values: median 0.537, IQR 0.368-0.710, range 0.014-0.936.
- NS2B median ρ 0.559; NS3 median ρ 0.527.
- Catalytic triad residue median ρ values: NS3:135=0.434, NS3:51=0.651, NS3:75=0.299.

## Figure 2 / Hierarchy

- Provisional gated-scale counts: residue 510, secondary_structure 399, chain 15; domain absent.
- Mean ρ rises from 0.527 at residue to 0.897 at secondary_structure; secondary_structure, motif and domain all have mean 0.897 and median 0.975.

## Figure 3 / Domains

- Catalytic Triad ρ 0.362-0.778, mean 0.526; coherence 0.572-0.965.
- Oxyanion Loop ρ 0.565-0.894, mean 0.680; coherence 0.712-0.889.
- NS3 unassigned ρ 0.903-0.978, mean 0.958.
- Domain mean-ρ ranking is in `numbers_audit.csv`.

## Figure 4 / Signed Effects

- Per-locus Figure 4 table: 924 loci tested; 302 signed; 152 FDR-significant.
- Catalytic Triad FDR-significant loci: 3; Oxyanion Loop: 11; combined: 14/152 (9.2%).
- Outside Catalytic Triad/Oxyanion Loop: 138/152 FDR-significant loci by count.

## Figure 5 / Cross-Serotype

- Conservation-class counts 0-4: 18, 49, 100, 63, 18.
- All-four-serotype reproducible positions: 18; catalytic-triad among them: 0.
- Median ρ by conservation class 0-4: 0.339, 0.395, 0.529, 0.625, 0.744.
- Exactly-three-serotype reproducible positions: 63.

## Figure 6 / Variance

- Median frac_tau2 0.925; IQR 0.833-0.963; range 0.358-0.994.
- Variance regimes: 28/32 replicate-dominated, 2/32 balanced, 2/32 sampling-dominated.
- Non-replicate-dominated exceptions: DENV1 120s Loop (sampling_dominated, frac_tau2=0.361), DENV1 Gly45 Turn (balanced, frac_tau2=0.594), DENV3 Gly45 Turn (balanced, frac_tau2=0.425), DENV4 B2b-C2 Hairpin (sampling_dominated, frac_tau2=0.358).
