# Claims To Edit

## Confirmed

- Per-serotype median residue-level ρ: 0.515 / 0.535 / 0.558 / 0.528 for DENV1-DENV4.
- n_loci range 199-247; n_gated_residue range 106-145.
- frac_mixed range 0.612-0.704.
- Catalytic Triad ρ 0.362-0.778, mean 0.526; Oxyanion Loop ρ 0.565-0.894, mean 0.680; NS3 unassigned ρ 0.903-0.978, mean 0.958.
- Catalytic Triad coherence 0.572-0.965; Oxyanion Loop coherence 0.712-0.889.
- Figure 5 conservation classes span 0-4 with counts 18, 49, 100, 63, 18; all-four core has 18 positions and 0 catalytic-triad residues.
- Figure 5C median ρ increases monotonically by conservation class.
- K=3, provisional ρ*=0.5, serotypes DENV1-DENV4.

## Wrong Or Needs Correction

- Variance wording: replace "29/32 replicate-dominated" with "28/32 replicate-dominated, with two balanced and two sampling-dominated entries" if using the stored `variance_regime` labels.
- Variance exceptions: do not call DENV3 Gly45 Turn sampling-dominated. It is `balanced` in `outputs_s4/variance_budget.parquet`; DENV1 Gly45 Turn is also balanced.
- Figure 5 caption: if the current rendered figure no longer marks catalytic-triad residues in Panel A, remove "with catalytic-triad residues marked."

## Denominator Clarification

- Figure 4/FDR claim should quote 152/924 FDR-significant loci from `outputs_s4/significance_screen.parquet`.
- The older 181 count comes from `outputs_s2/signed_screen.parquet` at ρ*=0.5 and reflects provisional screen pass / CI exclusion, not BH-FDR significance.
- "Catalytic minority" is supported by absolute Figure 4 counts: 14/152 FDR-significant loci are in Catalytic Triad or Oxyanion Loop. Do not imply a formal enrichment test.

## Too Vague / Softer Wording

- "Widespread significant effects" is supported, but state the denominator and avoid "enriched outside catalytic domains." More precise: "Most FDR-significant loci by count are outside the catalytic-triad and oxyanion-loop domains (138/152), although catalytic/oxyanion loci are not absent."
- "All but three entries are replicate-dominated" should become "Most entries are replicate-dominated (28/32); four are not, including two balanced and two sampling-dominated entries."
- Background claims about 65-79% identity and RMSD <1.5 A need external provenance; they are not verified by these frozen STRIDE outputs.
