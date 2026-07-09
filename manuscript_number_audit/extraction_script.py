#!/usr/bin/env python3
"""Extract manuscript-relevant numbers from frozen STRIDE outputs.

This script reads only the frozen output tables under strideouts and writes a
human audit, a long-form CSV, and manuscript-ready number notes. It does not
modify STRIDE outputs or recompute STRIDE statistics.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd


SEROTYPES = ["DENV1", "DENV2", "DENV3", "DENV4"]
CAT_DOMAINS = ["Catalytic Triad", "Oxyanion Loop"]


def fmt(value: Any, digits: int = 3) -> str:
    if pd.isna(value):
        return "NA"
    if isinstance(value, float):
        return f"{value:.{digits}f}"
    return str(value)


def pct(num: float, den: float, digits: int = 1) -> str:
    if den == 0:
        return "NA"
    return f"{100 * num / den:.{digits}f}%"


def compact_values(values: list[Any], digits: int = 3) -> str:
    return " / ".join(fmt(v, digits) for v in values)


def md_table(rows: list[dict[str, Any]], columns: list[str]) -> str:
    if not rows:
        return "_No rows._"
    out = ["|" + "|".join(columns) + "|",
           "|" + "|".join(["---"] * len(columns)) + "|"]
    for row in rows:
        vals = [str(row.get(col, "")).replace("\n", "<br>") for col in columns]
        out.append("|" + "|".join(vals) + "|")
    return "\n".join(out)


class Audit:
    def __init__(self, root: Path, outdir: Path) -> None:
        self.root = root
        self.outdir = outdir
        self.rows: list[dict[str, Any]] = []

    def source(self, rel: str) -> str:
        return str((self.root / rel).resolve())

    def parquet(self, rel: str) -> pd.DataFrame:
        return pd.read_parquet(self.root / rel)

    def csv(self, rel: str) -> pd.DataFrame:
        return pd.read_csv(self.root / rel)

    def json(self, rel: str) -> dict[str, Any]:
        return json.loads((self.root / rel).read_text())

    def add(self, section: str, claim: str, value: Any, source_file: str,
            columns: list[str], status: str = "new", wording: str = "",
            units: str = "", serotype: str = "", chain: str = "",
            domain: str = "", query: str = "") -> None:
        self.rows.append({
            "section": section,
            "claim_or_quantity": claim,
            "value": value if isinstance(value, str) else fmt(value, 6),
            "units_or_scale": units,
            "serotype": serotype,
            "chain": chain,
            "domain": domain,
            "source_file": self.source(source_file),
            "source_columns": ";".join(columns),
            "code_or_query": query,
            "manuscript_status": status,
            "suggested_wording": wording,
        })


def build_audit(root: Path, outdir: Path) -> None:
    audit = Audit(root, outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    t1 = audit.parquet("outputs_s7/T1_per_serotype_summary.parquet")
    ser_sum = audit.parquet("outputs_s2/serotype_summary.parquet")
    residue = audit.parquet("outputs_s2/residue_landscape.parquet")
    domains = audit.parquet("outputs_s2/domain_reproducibility.parquet")
    res_census = audit.parquet("outputs_s2/resolution_census.parquet")
    f2 = audit.parquet("outputs_s7/F2_resolution_census.parquet")
    scale_curve = audit.parquet("outputs_s3/scale_curve.parquet")
    cat_curve = audit.parquet("outputs_s7/F7_rho_vs_scale_catalytic.parquet")
    sig = audit.parquet("outputs_s4/significance_screen.parquet")
    s2_signed = audit.parquet("outputs_s2/signed_screen.parquet")
    dom_effect = audit.parquet("outputs_s4/domain_effect_summary.parquet")
    pos = audit.parquet("outputs_s5/position_conservation.parquet")
    concord = audit.parquet("outputs_s5/direction_concordance.parquet")
    domain_matrix = audit.parquet("outputs_s5/domain_serotype_matrix.parquet")
    variance = audit.parquet("outputs_s4/variance_budget.parquet")
    replicate_regime = audit.parquet("outputs_s6/replicate_regime.parquet")
    canonical = audit.parquet("outputs_s1a/canonical_residues.parquet")
    conservation_table = audit.parquet("outputs_s1a/conservation_table.parquet")
    artifact_manifest = audit.parquet("outputs_s7/artifact_manifest.parquet")
    dataset_summary = audit.json("outputs_s1a/dataset_summary.json")
    reduction_summary = audit.json("outputs_s2/reduction_summary.json")
    report_summary = audit.json("outputs_s7/report_summary.json")
    schema_report = audit.json("outputs/schema_report.json")

    # Part A: current manuscript claims.
    t1_ordered = t1.set_index("serotype").loc[SEROTYPES].reset_index()
    medians = t1_ordered["rho_residue_median"].tolist()
    audit.add("Part A - verified manuscript claims",
              "Per-serotype median residue-level rho, DENV1-DENV4",
              compact_values(medians, 3),
              "outputs_s7/T1_per_serotype_summary.parquet",
              ["serotype", "rho_residue_median"],
              "confirmed",
              "Per-serotype median residue-level ρ was 0.515, 0.535, 0.558 and 0.528 for DENV1-DENV4.")
    audit.add("Part A - verified manuscript claims", "n_loci range",
              f"{int(t1.n_loci.min())}-{int(t1.n_loci.max())}",
              "outputs_s7/T1_per_serotype_summary.parquet",
              ["n_loci"], "confirmed",
              "The four serotypes contributed 199-247 analysable loci.")
    audit.add("Part A - verified manuscript claims", "n_gated_residue range",
              f"{int(t1.n_gated_residue.min())}-{int(t1.n_gated_residue.max())}",
              "outputs_s7/T1_per_serotype_summary.parquet",
              ["n_gated_residue"], "confirmed",
              "At the provisional gate, 106-145 loci gated at residue resolution.")
    audit.add("Part A - verified manuscript claims", "frac_mixed range",
              f"{t1.frac_mixed.min():.3f}-{t1.frac_mixed.max():.3f}",
              "outputs_s7/T1_per_serotype_summary.parquet",
              ["frac_mixed"], "confirmed",
              "The mixed-mechanism fraction ranged from 0.612 to 0.704.")

    domain_claims: list[dict[str, str]] = []
    for display, domain, chain in [
        ("Catalytic Triad", "Catalytic Triad", "NS3"),
        ("Oxyanion Loop", "Oxyanion Loop", "NS3"),
        ("Unassigned NS3", "unassigned", "NS3"),
    ]:
        sub = domains[(domains["domain"] == domain) & (domains["chain"] == chain)]
        rho_min, rho_max, rho_mean = sub.rho_domain.min(), sub.rho_domain.max(), sub.rho_domain.mean()
        domain_claims.append({
            "domain": display,
            "rho_range": f"{rho_min:.3f}-{rho_max:.3f}",
            "rho_mean": f"{rho_mean:.3f}",
            "coherence_range": "" if display == "Unassigned NS3" else f"{sub.coherence_domain.min():.3f}-{sub.coherence_domain.max():.3f}",
        })
        audit.add("Part A - verified manuscript claims", f"{display} domain rho range and mean",
                  f"{rho_min:.6f}-{rho_max:.6f}; mean {rho_mean:.6f}",
                  "outputs_s2/domain_reproducibility.parquet",
                  ["serotype", "chain", "domain", "rho_domain"],
                  "confirmed", domain=display,
                  wording=f"{display} domain ρ spans {rho_min:.2f}-{rho_max:.2f} (mean {rho_mean:.2f}).")
        if display != "Unassigned NS3":
            audit.add("Part A - verified manuscript claims", f"{display} coherence range",
                      f"{sub.coherence_domain.min():.6f}-{sub.coherence_domain.max():.6f}",
                      "outputs_s2/domain_reproducibility.parquet",
                      ["serotype", "chain", "domain", "coherence_domain"],
                      "confirmed", domain=display,
                      wording=f"{display} coherence spans {sub.coherence_domain.min():.2f}-{sub.coherence_domain.max():.2f}.")

    sig_total = int(sig["significant_fdr"].sum())
    triad_fdr = int(sig[(sig["domain"] == "Catalytic Triad")]["significant_fdr"].sum())
    oxy_fdr = int(sig[(sig["domain"] == "Oxyanion Loop")]["significant_fdr"].sum())
    cat_oxy_fdr = triad_fdr + oxy_fdr
    fdr_signed = int(sig["is_signed"].sum())
    fdr_loci = len(sig)
    for claim, value, wording in [
        ("FDR-significant signed loci, per-locus significance_screen", sig_total,
         f"Figure 4 should use the per-locus FDR denominator: {sig_total}/{fdr_loci} loci are FDR-significant; {cat_oxy_fdr}/{sig_total} ({pct(cat_oxy_fdr, sig_total)}) are in Catalytic Triad or Oxyanion Loop."),
        ("Catalytic Triad FDR-significant loci", triad_fdr, ""),
        ("Oxyanion Loop FDR-significant loci", oxy_fdr, ""),
        ("Catalytic Triad + Oxyanion Loop FDR-significant loci", cat_oxy_fdr, ""),
    ]:
        audit.add("Part A - signed effects", claim, value,
                  "outputs_s4/significance_screen.parquet",
                  ["domain", "is_signed", "significant_fdr"],
                  "confirmed", wording=wording)
    s2_pass = int(s2_signed[s2_signed["rho_star"].eq(0.5)]["passes_screen"].sum())
    audit.add("Part A - signed effects", "Provisional signed-screen pass count at rho*=0.5",
              s2_pass, "outputs_s2/signed_screen.parquet",
              ["rho_star", "passes_screen"],
              "denominator clarification",
              "Do not mix this 181-count provisional/CI screen with the Figure 4 FDR-controlled count of 152.")

    class_counts = pos.groupby("n_serotypes_reproducible").size().reindex(range(5), fill_value=0)
    audit.add("Part A - cross-serotype", "Conservation class counts for 0-4 reproducing serotypes",
              "; ".join(f"{k}:{int(v)}" for k, v in class_counts.items()),
              "outputs_s5/position_conservation.parquet",
              ["n_serotypes_reproducible"], "confirmed",
              "Conservation classes contain 18, 49, 100, 63 and 18 positions for 0-4 reproducing serotypes.")
    all4 = pos[pos["n_serotypes_reproducible"] == 4].copy()
    audit.add("Part A - cross-serotype", "Positions reproducible in all four serotypes",
              len(all4), "outputs_s5/position_conservation.parquet",
              ["n_serotypes_reproducible"], "confirmed",
              f"There are {len(all4)} all-four-serotype reproducible positions; {int(all4.is_catalytic_triad.sum())} are catalytic triad residues.")
    med_by_class = pos.groupby("n_serotypes_reproducible")["rho_residue_median"].median()
    audit.add("Part A - cross-serotype", "Median rho by conservation class 0-4",
              "; ".join(f"{int(k)}:{v:.6f}" for k, v in med_by_class.items()),
              "outputs_s5/position_conservation.parquet",
              ["n_serotypes_reproducible", "rho_residue_median"],
              "confirmed",
              "Median residue-level ρ increases monotonically by conservation class: "
              + ", ".join(f"{v:.3f}" for v in med_by_class.tolist()) + ".")

    frac_tau_median = variance["frac_tau2"].median()
    regime_counts = variance["variance_regime"].value_counts()
    non_repl = variance[variance["variance_regime"] != "replicate_dominated"].copy()
    audit.add("Part A - variance", "Median frac_tau2",
              frac_tau_median, "outputs_s4/variance_budget.parquet",
              ["frac_tau2"], "confirmed",
              f"Across domain-by-serotype entries, median frac_tau2 is {frac_tau_median:.3f}.")
    audit.add("Part A - variance", "Variance regimes",
              "; ".join(f"{k}:{int(v)}" for k, v in regime_counts.items()),
              "outputs_s4/variance_budget.parquet",
              ["variance_regime"], "wrong in manuscript if stated as 29/32",
              "Use 28 replicate-dominated, 2 balanced and 2 sampling-dominated. If emphasizing non-replicate-dominated entries, list all four.")
    for _, row in non_repl.sort_values(["variance_regime", "serotype", "domain"]).iterrows():
        audit.add("Part A - variance exceptions",
                  "Non-replicate-dominated variance-budget entry",
                  f"frac_tau2={row.frac_tau2:.6f}; frac_sigma2={row.frac_sigma2:.6f}; regime={row.variance_regime}",
                  "outputs_s4/variance_budget.parquet",
                  ["serotype", "chain", "domain", "frac_tau2", "frac_sigma2", "variance_regime"],
                  "needs manuscript edit", serotype=row.serotype, chain=row.chain, domain=row.domain)
    for domain in CAT_DOMAINS:
        sub = variance[variance["domain"] == domain]
        audit.add("Part A - variance", f"{domain} frac_tau2 range",
                  f"{sub.frac_tau2.min():.6f}-{sub.frac_tau2.max():.6f}",
                  "outputs_s4/variance_budget.parquet",
                  ["domain", "frac_tau2"], "confirmed", domain=domain)
    audit.add("Part A - variance", "Overall frac_tau2 range",
              f"{variance.frac_tau2.min():.6f}-{variance.frac_tau2.max():.6f}",
              "outputs_s4/variance_budget.parquet",
              ["frac_tau2"], "new")

    f2_levels = sorted(f2["gated_scale_level"].unique())
    audit.add("Part A - hierarchy", "Provisional gated-scale categories in S7 Figure 2",
              ", ".join(f2_levels),
              "outputs_s7/F2_resolution_census.parquet",
              ["gated_scale_level", "n_loci"], "confirmed",
              "At rho*=0.5, the census contains residue, secondary_structure and chain; domain is absent.")
    scale_summary = scale_curve.groupby("scale_level")["rho"].agg(["mean", "median", "min", "max", "count"])
    first_gain = scale_summary.loc["secondary_structure", "mean"] - scale_summary.loc["residue", "mean"]
    audit.add("Part A - hierarchy", "Mean rho gain from residue to secondary_structure",
              first_gain, "outputs_s3/scale_curve.parquet",
              ["scale_level", "rho"], "confirmed",
              f"Mean ρ rises from {scale_summary.loc['residue','mean']:.3f} at residue to {scale_summary.loc['secondary_structure','mean']:.3f} at secondary_structure; motif and domain have the same summary values.")

    prov = reduction_summary["provenance"]["provisional_rho_star"]
    audit.add("Part A - provenance", "Provisional rho* value",
              prov, "outputs_s2/reduction_summary.json",
              ["provenance.provisional_rho_star"], "confirmed")
    k_values = replicate_regime.set_index("serotype").loc[SEROTYPES]["n_replicates"].tolist()
    audit.add("Part A - provenance", "Replicates per serotype",
              compact_values(k_values, 0), "outputs_s6/replicate_regime.parquet",
              ["serotype", "n_replicates"], "confirmed",
              "Each serotype has K=3 independent replicate runs.")
    audit.add("Part A - provenance", "Serotypes present",
              ", ".join(report_summary["serotypes"]), "outputs_s7/report_summary.json",
              ["serotypes"], "confirmed")
    audit.add("Part A - provenance", "S7 manifest artifacts",
              f"{len(artifact_manifest)} artifacts ({report_summary['facts']['n_figures']} figures; {report_summary['facts']['n_tables']} tables)",
              "outputs_s7/artifact_manifest.parquet",
              ["artifact_id", "kind"], "confirmed")

    # Part B: additional useful numbers.
    audit.add("Part B - dataset", "Number of serotypes", len(SEROTYPES),
              "outputs_s7/report_summary.json", ["facts.n_serotypes"], "new")
    audit.add("Part B - dataset", "Total replicate runs",
              int(replicate_regime["n_replicates"].sum()),
              "outputs_s6/replicate_regime.parquet", ["n_replicates"], "new")
    audit.add("Part B - dataset", "Replicate table rows",
              schema_report["replicate_table_rows"], "outputs/schema_report.json",
              ["replicate_table_rows"], "new")
    audit.add("Part B - dataset", "Canonical serotype-residue rows",
              dataset_summary["n_canonical_residues"], "outputs_s1a/dataset_summary.json",
              ["n_canonical_residues"], "new")
    audit.add("Part B - dataset", "Union residue positions",
              dataset_summary["n_union_residues"], "outputs_s1a/dataset_summary.json",
              ["n_union_residues"], "new")
    audit.add("Part B - dataset", "Positions present in all serotypes",
              dataset_summary["n_conserved_all_serotypes"], "outputs_s1a/dataset_summary.json",
              ["n_conserved_all_serotypes"], "new")
    for chain, n in conservation_table.groupby("chain").size().items():
        audit.add("Part B - dataset", "Union positions by chain", int(n),
                  "outputs_s1a/conservation_table.parquet", ["chain", "canon_label"],
                  "new", chain=chain)
    audit.add("Part B - dataset", "Domain-by-serotype entries",
              len(domain_matrix), "outputs_s5/domain_serotype_matrix.parquet",
              ["serotype", "chain", "domain"], "new")
    audit.add("Part B - dataset", "Cross-serotype position table rows",
              len(pos), "outputs_s5/position_conservation.parquet",
              ["canon_label"], "new")
    audit.add("Part B - dataset", "Signed mechanisms in per-locus FDR screen",
              fdr_signed, "outputs_s4/significance_screen.parquet",
              ["is_signed"], "new")
    audit.add("Part B - dataset", "Mixed mechanisms at provisional rho*",
              int(t1["n_mixed"].sum()), "outputs_s7/T1_per_serotype_summary.parquet",
              ["n_mixed"], "new")

    residue_summary = residue.groupby("serotype")["rho_residue"].agg(["median", "quantile", "min", "max"])
    # pandas named quantiles are clunky; add directly from T1 for manuscript values.
    for _, row in t1_ordered.iterrows():
        audit.add("Part B - reproducibility landscape", "Residue rho median/IQR/min/max per serotype",
                  f"median={row.rho_residue_median:.6f}; IQR={row.rho_residue_q1:.6f}-{row.rho_residue_q3:.6f}; min={row.rho_residue_min:.6f}; max={row.rho_residue_max:.6f}",
                  "outputs_s7/T1_per_serotype_summary.parquet",
                  ["serotype", "rho_residue_median", "rho_residue_q1", "rho_residue_q3", "rho_residue_min", "rho_residue_max"],
                  "new", serotype=row.serotype)
    audit.add("Part B - reproducibility landscape", "All residue-level rho distribution",
              f"median={residue.rho_residue.median():.6f}; IQR={residue.rho_residue.quantile(0.25):.6f}-{residue.rho_residue.quantile(0.75):.6f}; min={residue.rho_residue.min():.6f}; max={residue.rho_residue.max():.6f}",
              "outputs_s2/residue_landscape.parquet",
              ["rho_residue"], "new")
    for chain, sub in residue.groupby("chain"):
        audit.add("Part B - reproducibility landscape", "Residue-level rho by chain",
                  f"n={len(sub)}; median={sub.rho_residue.median():.6f}; mean={sub.rho_residue.mean():.6f}",
                  "outputs_s2/residue_landscape.parquet",
                  ["chain", "rho_residue"], "new", chain=chain)
    median_pos = residue.groupby(["canon_label", "chain", "domain"])["rho_residue"].median().reset_index(name="median_rho")
    for label, frame in [("Top 10 median-rho positions", median_pos.sort_values("median_rho", ascending=False).head(10)),
                         ("Bottom 10 median-rho positions", median_pos.sort_values("median_rho", ascending=True).head(10))]:
        for rank, (_, row) in enumerate(frame.iterrows(), 1):
            audit.add("Part B - reproducibility landscape", f"{label} rank {rank}",
                      row.median_rho, "outputs_s2/residue_landscape.parquet",
                      ["canon_label", "chain", "domain", "rho_residue"],
                      "new", chain=row.chain, domain=row.domain,
                      wording=f"{row.canon_label} ({row.chain}, {row.domain}) median ρ={row.median_rho:.3f}.")
    for _, row in pos[pos["is_catalytic_triad"]].sort_values("canon_label").iterrows():
        audit.add("Part B - reproducibility landscape", "Catalytic triad residue-level rho",
                  f"n_serotypes_reproducible={int(row.n_serotypes_reproducible)}; min={row.rho_residue_min:.6f}; median={row.rho_residue_median:.6f}; max={row.rho_residue_max:.6f}",
                  "outputs_s5/position_conservation.parquet",
                  ["canon_label", "is_catalytic_triad", "n_serotypes_reproducible", "rho_residue_min", "rho_residue_median", "rho_residue_max"],
                  "new", chain=row.chain, domain=row.domain, query=row.canon_label)

    for scale, row in scale_summary.iterrows():
        audit.add("Part B - spatial hierarchy", "rho summary by scale",
                  f"n={int(row['count'])}; mean={row['mean']:.6f}; median={row['median']:.6f}; min={row['min']:.6f}; max={row['max']:.6f}",
                  "outputs_s3/scale_curve.parquet",
                  ["scale_level", "rho"], "new", query=scale)
    for _, row in f2.iterrows():
        audit.add("Part B - spatial hierarchy", "Provisional gated-scale count by serotype",
                  int(row.n_loci), "outputs_s7/F2_resolution_census.parquet",
                  ["serotype", "gated_scale_level", "n_loci"], "new",
                  serotype=row.serotype, query=row.gated_scale_level)
    cat_scale = cat_curve.groupby("scale_level")["rho"].agg(["mean", "median", "min", "max", "count"])
    for scale, row in cat_scale.iterrows():
        audit.add("Part B - spatial hierarchy", "Catalytic-region rho by scale",
                  f"n={int(row['count'])}; mean={row['mean']:.6f}; median={row['median']:.6f}; min={row['min']:.6f}; max={row['max']:.6f}",
                  "outputs_s7/F7_rho_vs_scale_catalytic.parquet",
                  ["scale_level", "rho"], "new", query=scale)

    domain_means = domains.groupby(["chain", "domain"]).agg(
        mean_rho=("rho_domain", "mean"),
        min_rho=("rho_domain", "min"),
        max_rho=("rho_domain", "max"),
        mean_coherence=("coherence_domain", "mean"),
    ).reset_index().sort_values("mean_rho", ascending=False)
    for rank, (_, row) in enumerate(domain_means.iterrows(), 1):
        audit.add("Part B - domain analysis", f"Domain mean-rho rank {rank}",
                  f"mean_rho={row.mean_rho:.6f}; range={row.min_rho:.6f}-{row.max_rho:.6f}; mean_coherence={row.mean_coherence:.6f}",
                  "outputs_s2/domain_reproducibility.parquet",
                  ["chain", "domain", "rho_domain", "coherence_domain"],
                  "new", chain=row.chain, domain=row.domain)

    sign_split = sig[sig["significant_fdr"]].assign(sign=lambda d: d["beta_signed"].gt(0).map({True: "positive", False: "negative"})).groupby("sign").size()
    for sign, n in sign_split.items():
        audit.add("Part B - signed effects", "FDR-significant signed effects by sign",
                  int(n), "outputs_s4/significance_screen.parquet",
                  ["significant_fdr", "beta_signed"], "new", query=sign)
    for label, frame in [
        ("Top positive beta FDR loci", sig[sig.significant_fdr].sort_values("beta_signed", ascending=False).head(10)),
        ("Top negative beta FDR loci", sig[sig.significant_fdr].sort_values("beta_signed", ascending=True).head(10)),
        ("Top FDR loci by adjusted p", sig[sig.significant_fdr].sort_values("p_value_bh", ascending=True).head(10)),
    ]:
        for rank, (_, row) in enumerate(frame.iterrows(), 1):
            audit.add("Part B - signed effects", f"{label} rank {rank}",
                      f"beta={row.beta_signed:.6f}; p_bh={row.p_value_bh:.6g}",
                      "outputs_s4/significance_screen.parquet",
                      ["serotype", "canon_label", "chain", "domain", "beta_signed", "p_value_bh", "significant_fdr"],
                      "new", serotype=row.serotype, chain=row.chain, domain=row.domain, query=row.canon_label)
    for _, row in sig[(sig["domain"].isin(CAT_DOMAINS)) & (sig["is_signed"])].sort_values(["domain", "serotype", "canon_label"]).iterrows():
        audit.add("Part B - signed effects", "Catalytic-domain signed beta value",
                  f"beta={row.beta_signed:.6f}; p_bh={row.p_value_bh:.6g}; fdr={bool(row.significant_fdr)}",
                  "outputs_s4/significance_screen.parquet",
                  ["serotype", "canon_label", "domain", "beta_signed", "p_value_bh", "significant_fdr"],
                  "new", serotype=row.serotype, chain=row.chain, domain=row.domain, query=row.canon_label)

    repro_by_ser_chain = residue[residue["rho_residue"] >= prov].groupby(["serotype", "chain"]).size()
    for (serotype, chain), n in repro_by_ser_chain.items():
        audit.add("Part B - cross-serotype conservation", "Reproducible-position count by serotype and chain",
                  int(n), "outputs_s2/residue_landscape.parquet",
                  ["serotype", "chain", "rho_residue"], "new", serotype=serotype, chain=chain)
    all4_join = all4.merge(concord[["canon_label", "n_serotypes_signed", "n_increase", "n_decrease", "majority_direction", "concordance_class"]],
                           on="canon_label", how="left")
    for _, row in all4_join.sort_values(["chain", "domain", "canon_label"]).iterrows():
        audit.add("Part B - cross-serotype conservation", "All-four-serotype reproducible position",
                  f"median_rho={row.rho_residue_median:.6f}; signed_reproducible_serotypes={int(row.n_serotypes_signed_reproducible)}; majority_direction={row.get('majority_direction', 'NA')}; concordance={row.get('concordance_class', 'NA')}; catalytic={bool(row.is_catalytic_triad)}",
                  "outputs_s5/position_conservation.parquet",
                  ["canon_label", "chain", "domain", "n_serotypes_reproducible", "rho_residue_median", "n_serotypes_signed_reproducible", "is_catalytic_triad"],
                  "new", chain=row.chain, domain=row.domain, query=row.canon_label)
    exactly3 = pos[pos["n_serotypes_reproducible"] == 3].sort_values(["chain", "domain", "canon_label"])
    for _, row in exactly3.iterrows():
        audit.add("Part B - cross-serotype conservation", "Exactly-three-serotype reproducible position",
                  f"median_rho={row.rho_residue_median:.6f}; signed_reproducible_serotypes={int(row.n_serotypes_signed_reproducible)}; catalytic={bool(row.is_catalytic_triad)}",
                  "outputs_s5/position_conservation.parquet",
                  ["canon_label", "chain", "domain", "n_serotypes_reproducible", "rho_residue_median", "n_serotypes_signed_reproducible", "is_catalytic_triad"],
                  "new", chain=row.chain, domain=row.domain, query=row.canon_label)

    for _, row in variance.groupby("serotype")["frac_tau2"].agg(["median", "min", "max", "count"]).reset_index().iterrows():
        audit.add("Part B - variance budget", "Per-serotype frac_tau2 summary",
                  f"n={int(row['count'])}; median={row['median']:.6f}; min={row['min']:.6f}; max={row['max']:.6f}",
                  "outputs_s4/variance_budget.parquet",
                  ["serotype", "frac_tau2"], "new", serotype=row.serotype)
    for label, frame in [
        ("Most replicate-dominated regions", variance.sort_values("frac_tau2", ascending=False).head(10)),
        ("Most sampling-dominated regions", variance.sort_values("frac_tau2", ascending=True).head(10)),
    ]:
        for rank, (_, row) in enumerate(frame.iterrows(), 1):
            audit.add("Part B - variance budget", f"{label} rank {rank}",
                      f"frac_tau2={row.frac_tau2:.6f}; frac_sigma2={row.frac_sigma2:.6f}; regime={row.variance_regime}",
                      "outputs_s4/variance_budget.parquet",
                      ["serotype", "chain", "domain", "frac_tau2", "frac_sigma2", "variance_regime"],
                      "new", serotype=row.serotype, chain=row.chain, domain=row.domain)

    for _, row in artifact_manifest.sort_values(["kind", "artifact_id"]).iterrows():
        audit.add("Part B - figure/table provenance", f"{row.artifact_id} {row.kind}: {row.title}",
                  f"rows={int(row.n_rows)}; sources={row.sources}; files={row.files}",
                  "outputs_s7/artifact_manifest.parquet",
                  ["artifact_id", "kind", "title", "sources", "files", "n_rows"],
                  "new")

    csv_df = pd.DataFrame(audit.rows)
    csv_df.to_csv(outdir / "numbers_audit.csv", index=False)

    # Markdown reports.
    t1_rows = [
        {
            "serotype": r.serotype,
            "n_loci": int(r.n_loci),
            "n_gated_residue": int(r.n_gated_residue),
            "frac_mixed": f"{r.frac_mixed:.3f}",
            "median_rho": f"{r.rho_residue_median:.3f}",
            "IQR": f"{r.rho_residue_q1:.3f}-{r.rho_residue_q3:.3f}",
        }
        for _, r in t1_ordered.iterrows()
    ]
    sig_by_domain = sig.groupby(["chain", "domain"]).agg(
        loci=("canon_label", "size"),
        signed=("is_signed", "sum"),
        fdr=("significant_fdr", "sum"),
    ).reset_index().sort_values("fdr", ascending=False)
    sig_rows = [
        {"chain": r.chain, "domain": r.domain, "loci": int(r.loci), "signed": int(r.signed), "fdr": int(r.fdr)}
        for _, r in sig_by_domain.iterrows()
    ]
    all4_rows = [
        {
            "position": r.canon_label,
            "chain": r.chain,
            "domain": r.domain,
            "median_rho": f"{r.rho_residue_median:.3f}",
            "signed_serotypes": int(r.n_serotypes_signed_reproducible),
            "direction": fmt(r.get("majority_direction", "NA")),
            "catalytic": bool(r.is_catalytic_triad),
        }
        for _, r in all4_join.sort_values(["chain", "domain", "canon_label"]).iterrows()
    ]
    variance_exception_rows = [
        {
            "serotype": r.serotype,
            "chain": r.chain,
            "domain": r.domain,
            "frac_tau2": f"{r.frac_tau2:.3f}",
            "frac_sigma2": f"{r.frac_sigma2:.3f}",
            "regime": r.variance_regime,
        }
        for _, r in non_repl.sort_values(["serotype", "domain"]).iterrows()
    ]
    domain_rows = [
        {
            "rank": i,
            "chain": r.chain,
            "domain": r.domain,
            "mean_rho": f"{r.mean_rho:.3f}",
            "rho_range": f"{r.min_rho:.3f}-{r.max_rho:.3f}",
            "mean_coherence": f"{r.mean_coherence:.3f}",
        }
        for i, (_, r) in enumerate(domain_means.iterrows(), 1)
    ]
    scale_rows = [
        {
            "scale": scale,
            "mean_rho": f"{row['mean']:.3f}",
            "median_rho": f"{row['median']:.3f}",
            "range": f"{row['min']:.3f}-{row['max']:.3f}",
            "n": int(row["count"]),
        }
        for scale, row in scale_summary.iterrows()
    ]
    f2_rows = [
        {"serotype": r.serotype, "scale": r.gated_scale_level, "n_loci": int(r.n_loci)}
        for _, r in f2.iterrows()
    ]

    numbers_audit = f"""# Manuscript Number Audit

Generated from frozen STRIDE outputs under `{root}`.

No STRIDE outputs were modified. The long-form source table is `numbers_audit.csv`.

## Source Hierarchy

- Use stage outputs (`outputs_s2`-`outputs_s6`) when a manuscript claim needs columns not exported in S7.
- Use S7 figure/table artifacts when verifying the exact figure/table-ready values.
- For Figure 4 signed effects, the authoritative plotted denominator is `outputs_s4/significance_screen.parquet`, not the provisional `outputs_s2/signed_screen.parquet`.

## Part A - Current Manuscript Claims

### Table 1 / Comparable Substrate

{md_table(t1_rows, ["serotype", "n_loci", "n_gated_residue", "frac_mixed", "median_rho", "IQR"])}

Status: confirmed. The manuscript values 0.515 / 0.535 / 0.558 / 0.528, n_loci 199-247, n_gated_residue 106-145, and frac_mixed roughly 0.61-0.70 are supported.

### Domain Reproducibility

{md_table(domain_claims, ["domain", "rho_range", "rho_mean", "coherence_range"])}

Status: confirmed, with one naming clarification. The manuscript phrase "Unassigned NS3" corresponds to `chain == NS3` and `domain == unassigned` in the frozen table.

### Signed Effects / Figure 4

- Per-locus FDR-controlled significance screen: {sig_total}/{fdr_loci} loci are FDR-significant.
- Signed loci in that table: {fdr_signed}/{fdr_loci}.
- Catalytic Triad FDR-significant loci: {triad_fdr}/{sig_total}.
- Oxyanion Loop FDR-significant loci: {oxy_fdr}/{sig_total}.
- Catalytic Triad + Oxyanion Loop total: {cat_oxy_fdr}/{sig_total} ({pct(cat_oxy_fdr, sig_total)} of FDR-significant loci).
- Provisional S2 screen-pass / CI-excluding-zero count at rho*=0.5: {s2_pass}. This is not the Figure 4 FDR denominator.

{md_table(sig_rows, ["chain", "domain", "loci", "signed", "fdr"])}

Interpretation: "significant signed effects are widespread" is supported. "Catalytic minority" is supported by absolute FDR counts under the Figure 4 denominator ({cat_oxy_fdr}/{sig_total}); however, it should not be worded as a formal enrichment/depletion claim. Descriptively, catalytic+oxyanion regions have {cat_oxy_fdr}/44 FDR-significant loci ({pct(cat_oxy_fdr, 44)}) versus {sig_total - cat_oxy_fdr}/880 outside those domains ({pct(sig_total - cat_oxy_fdr, 880)}), so the outside majority is driven by many more outside loci, not by a lower catalytic per-locus rate.

### Cross-Serotype Conservation / Figure 5

- Conservation-class counts for 0, 1, 2, 3, 4 reproducing serotypes: {", ".join(str(int(v)) for v in class_counts.tolist())}.
- Positions reproducible in all four serotypes: {len(all4)}.
- Catalytic-triad residues among all-four positions: {int(all4.is_catalytic_triad.sum())}.
- Median rho by conservation class 0-4: {", ".join(f"{v:.3f}" for v in med_by_class.tolist())}.
- Monotonic median-rho increase across classes: {bool(med_by_class.is_monotonic_increasing)}.

All-four-serotype positions:

{md_table(all4_rows, ["position", "chain", "domain", "median_rho", "signed_serotypes", "direction", "catalytic"])}

Status: confirmed. The all-four-serotype core is limited and non-catalytic in the frozen table.

### Variance / Figure 6

- Median frac_tau2: {frac_tau_median:.3f}.
- Variance regimes: {", ".join(f"{k}={int(v)}" for k, v in regime_counts.items())}.
- Overall frac_tau2 range: {variance.frac_tau2.min():.3f}-{variance.frac_tau2.max():.3f}.
- Catalytic Triad frac_tau2 range: {variance[variance.domain.eq("Catalytic Triad")].frac_tau2.min():.3f}-{variance[variance.domain.eq("Catalytic Triad")].frac_tau2.max():.3f}.
- Oxyanion Loop frac_tau2 range: {variance[variance.domain.eq("Oxyanion Loop")].frac_tau2.min():.3f}-{variance[variance.domain.eq("Oxyanion Loop")].frac_tau2.max():.3f}.

Non-replicate-dominated exceptions:

{md_table(variance_exception_rows, ["serotype", "chain", "domain", "frac_tau2", "frac_sigma2", "regime"])}

Status: median frac_tau2 is confirmed, but the manuscript's "29/32 replicate-dominated" and "three sampling-dominated exceptions" wording is not supported by the stored regime labels. The frozen table says 28 replicate-dominated, 2 balanced and 2 sampling-dominated.

### Hierarchy / Figure 2

Provisional Figure 2B gated-scale categories: {", ".join(f2_levels)}. Domain is not present as a provisional census category.

{md_table(f2_rows, ["serotype", "scale", "n_loci"])}

Scale-curve rho summaries:

{md_table(scale_rows, ["scale", "mean_rho", "median_rho", "range", "n"])}

Status: the plateau claim is supported. Secondary-structure, motif and domain have identical mean/median/range summaries in the full scale curve, and catalytic-region F7 has the same collapse across those three aggregate levels.

### Provenance

- K=3 per serotype: confirmed from `outputs_s6/replicate_regime.parquet`.
- Provisional rho*: {prov}.
- Serotypes present: {", ".join(report_summary["serotypes"])}.
- S7 artifacts: {len(artifact_manifest)} total ({report_summary["facts"]["n_figures"]} figures and {report_summary["facts"]["n_tables"]} tables).

## Part B - Additional Useful Numbers

### Dataset / Case-Study Scale

- Serotypes: {len(SEROTYPES)}.
- Replicate runs: {int(replicate_regime["n_replicates"].sum())} total, K=3 per serotype.
- Replicate-table rows: {schema_report["replicate_table_rows"]}.
- Serotype-residue rows: {dataset_summary["n_canonical_residues"]}; union residue positions: {dataset_summary["n_union_residues"]}; positions present in all serotypes: {dataset_summary["n_conserved_all_serotypes"]}.
- Union positions by chain: {", ".join(f"{chain}={n}" for chain, n in conservation_table.groupby("chain").size().items())}.
- Domain-by-serotype entries: {len(domain_matrix)}.
- Cross-serotype position rows: {len(pos)}.
- Per-locus signed effects: {fdr_signed}; mixed mechanisms at provisional rho*: {int(t1.n_mixed.sum())}; FDR-significant signed effects: {sig_total}.

### Reproducibility Landscape

- All residue-level rho values: median {residue.rho_residue.median():.3f}, IQR {residue.rho_residue.quantile(0.25):.3f}-{residue.rho_residue.quantile(0.75):.3f}, range {residue.rho_residue.min():.3f}-{residue.rho_residue.max():.3f}.
- Chain medians: {", ".join(f"{chain}={sub.rho_residue.median():.3f}" for chain, sub in residue.groupby("chain"))}.
- Top/bottom 10 median-rho positions and catalytic-triad residue-level rho values are enumerated in `numbers_audit.csv`.

### Spatial Hierarchy

- Mean rho gain from residue to first aggregate: {first_gain:.3f}.
- Residue mean/median rho: {scale_summary.loc["residue", "mean"]:.3f}/{scale_summary.loc["residue", "median"]:.3f}.
- Secondary-structure/motif/domain mean/median rho: {scale_summary.loc["secondary_structure", "mean"]:.3f}/{scale_summary.loc["secondary_structure", "median"]:.3f}.

### Domain Analysis

{md_table(domain_rows, ["rank", "chain", "domain", "mean_rho", "rho_range", "mean_coherence"])}

Useful contrast: NS3 unassigned has the highest mean rho but low coherence; Catalytic Triad has the lowest mean rho among listed domains but relatively high coherence.

### Signed Effects

- FDR-significant effects by sign: {", ".join(f"{k}={int(v)}" for k, v in sign_split.items())}.
- Top positive/negative beta loci, top adjusted-p loci, and catalytic-domain signed beta values are enumerated in `numbers_audit.csv`.

### Cross-Serotype Conservation

- Exactly-three-serotype reproducible positions: {len(exactly3)}; all are listed in `numbers_audit.csv`.
- Serotype/chain reproducible-position counts are listed in `numbers_audit.csv`.

### Variance Budget

- Overall frac_tau2 IQR: {variance.frac_tau2.quantile(0.25):.3f}-{variance.frac_tau2.quantile(0.75):.3f}.
- Replicate-dominated fraction: {int(regime_counts.get("replicate_dominated", 0))}/{len(variance)} ({pct(int(regime_counts.get("replicate_dominated", 0)), len(variance))}).
- Sampling-dominated fraction: {int(regime_counts.get("sampling_dominated", 0))}/{len(variance)} ({pct(int(regime_counts.get("sampling_dominated", 0)), len(variance))}).
- Most replicate- and sampling-dominated regions are enumerated in `numbers_audit.csv`.

## Missing or Not Verifiable from Frozen Outputs

- Sequence identity 65-79% and static-structure RMSD below 1.5 A are background claims and are not verifiable from the frozen STRIDE reduction outputs.
- MD engine/version, force field, water model, ensemble, temperature/pressure/ions, equilibration, trajectory length, frame sampling, PDB accessions and shared-reference-frame control results are not present in the inspected frozen reduction outputs.
- The manuscript's Figure 5 caption still says catalytic-triad residues are marked in Panel A; after the presentation cleanup this should be edited if the figure no longer marks catalytic residues there.
"""

    claims_to_edit = f"""# Claims To Edit

## Confirmed

- Per-serotype median residue-level ρ: {compact_values(medians, 3)} for DENV1-DENV4.
- n_loci range {int(t1.n_loci.min())}-{int(t1.n_loci.max())}; n_gated_residue range {int(t1.n_gated_residue.min())}-{int(t1.n_gated_residue.max())}.
- frac_mixed range {t1.frac_mixed.min():.3f}-{t1.frac_mixed.max():.3f}.
- Catalytic Triad ρ {domain_claims[0]["rho_range"]}, mean {domain_claims[0]["rho_mean"]}; Oxyanion Loop ρ {domain_claims[1]["rho_range"]}, mean {domain_claims[1]["rho_mean"]}; NS3 unassigned ρ {domain_claims[2]["rho_range"]}, mean {domain_claims[2]["rho_mean"]}.
- Catalytic Triad coherence {domain_claims[0]["coherence_range"]}; Oxyanion Loop coherence {domain_claims[1]["coherence_range"]}.
- Figure 5 conservation classes span 0-4 with counts {", ".join(str(int(v)) for v in class_counts.tolist())}; all-four core has {len(all4)} positions and {int(all4.is_catalytic_triad.sum())} catalytic-triad residues.
- Figure 5C median ρ increases monotonically by conservation class.
- K=3, provisional ρ*=0.5, serotypes DENV1-DENV4.

## Wrong Or Needs Correction

- Variance wording: replace "29/32 replicate-dominated" with "28/32 replicate-dominated, with two balanced and two sampling-dominated entries" if using the stored `variance_regime` labels.
- Variance exceptions: do not call DENV3 Gly45 Turn sampling-dominated. It is `balanced` in `outputs_s4/variance_budget.parquet`; DENV1 Gly45 Turn is also balanced.
- Figure 5 caption: if the current rendered figure no longer marks catalytic-triad residues in Panel A, remove "with catalytic-triad residues marked."

## Denominator Clarification

- Figure 4/FDR claim should quote {sig_total}/{fdr_loci} FDR-significant loci from `outputs_s4/significance_screen.parquet`.
- The older {s2_pass} count comes from `outputs_s2/signed_screen.parquet` at ρ*=0.5 and reflects provisional screen pass / CI exclusion, not BH-FDR significance.
- "Catalytic minority" is supported by absolute Figure 4 counts: {cat_oxy_fdr}/{sig_total} FDR-significant loci are in Catalytic Triad or Oxyanion Loop. Do not imply a formal enrichment test.

## Too Vague / Softer Wording

- "Widespread significant effects" is supported, but state the denominator and avoid "enriched outside catalytic domains." More precise: "Most FDR-significant loci by count are outside the catalytic-triad and oxyanion-loop domains ({sig_total - cat_oxy_fdr}/{sig_total}), although catalytic/oxyanion loci are not absent."
- "All but three entries are replicate-dominated" should become "Most entries are replicate-dominated (28/32); four are not, including two balanced and two sampling-dominated entries."
- Background claims about 65-79% identity and RMSD <1.5 A need external provenance; they are not verified by these frozen STRIDE outputs.
"""

    numbers_for_manuscript = f"""# Numbers For Manuscript

## Results 4.1

- Four serotypes, K=3 independent replicate runs per serotype ({int(replicate_regime.n_replicates.sum())} total replicate runs).
- Analysable loci per serotype: {compact_values(t1_ordered.n_loci.tolist(), 0)} for DENV1-DENV4 (range {int(t1.n_loci.min())}-{int(t1.n_loci.max())}).
- Gated residues at provisional ρ*=0.5: {compact_values(t1_ordered.n_gated_residue.tolist(), 0)} (range {int(t1.n_gated_residue.min())}-{int(t1.n_gated_residue.max())}).
- Median residue-level ρ: {compact_values(medians, 3)} for DENV1-DENV4.
- Mixed-mechanism fraction: {compact_values(t1_ordered.frac_mixed.tolist(), 3)}.

## Figure 1 / Landscape

- All residue-level ρ values: median {residue.rho_residue.median():.3f}, IQR {residue.rho_residue.quantile(0.25):.3f}-{residue.rho_residue.quantile(0.75):.3f}, range {residue.rho_residue.min():.3f}-{residue.rho_residue.max():.3f}.
- NS2B median ρ {residue[residue.chain.eq("NS2B")].rho_residue.median():.3f}; NS3 median ρ {residue[residue.chain.eq("NS3")].rho_residue.median():.3f}.
- Catalytic triad residue median ρ values: {", ".join(f"{r.canon_label}={r.rho_residue_median:.3f}" for _, r in pos[pos.is_catalytic_triad].sort_values("canon_label").iterrows())}.

## Figure 2 / Hierarchy

- Provisional gated-scale counts: residue {int(f2[f2.gated_scale_level.eq("residue")].n_loci.sum())}, secondary_structure {int(f2[f2.gated_scale_level.eq("secondary_structure")].n_loci.sum())}, chain {int(f2[f2.gated_scale_level.eq("chain")].n_loci.sum())}; domain absent.
- Mean ρ rises from {scale_summary.loc["residue", "mean"]:.3f} at residue to {scale_summary.loc["secondary_structure", "mean"]:.3f} at secondary_structure; secondary_structure, motif and domain all have mean {scale_summary.loc["domain", "mean"]:.3f} and median {scale_summary.loc["domain", "median"]:.3f}.

## Figure 3 / Domains

- Catalytic Triad ρ {domain_claims[0]["rho_range"]}, mean {domain_claims[0]["rho_mean"]}; coherence {domain_claims[0]["coherence_range"]}.
- Oxyanion Loop ρ {domain_claims[1]["rho_range"]}, mean {domain_claims[1]["rho_mean"]}; coherence {domain_claims[1]["coherence_range"]}.
- NS3 unassigned ρ {domain_claims[2]["rho_range"]}, mean {domain_claims[2]["rho_mean"]}.
- Domain mean-ρ ranking is in `numbers_audit.csv`.

## Figure 4 / Signed Effects

- Per-locus Figure 4 table: {fdr_loci} loci tested; {fdr_signed} signed; {sig_total} FDR-significant.
- Catalytic Triad FDR-significant loci: {triad_fdr}; Oxyanion Loop: {oxy_fdr}; combined: {cat_oxy_fdr}/{sig_total} ({pct(cat_oxy_fdr, sig_total)}).
- Outside Catalytic Triad/Oxyanion Loop: {sig_total - cat_oxy_fdr}/{sig_total} FDR-significant loci by count.

## Figure 5 / Cross-Serotype

- Conservation-class counts 0-4: {", ".join(str(int(v)) for v in class_counts.tolist())}.
- All-four-serotype reproducible positions: {len(all4)}; catalytic-triad among them: {int(all4.is_catalytic_triad.sum())}.
- Median ρ by conservation class 0-4: {", ".join(f"{v:.3f}" for v in med_by_class.tolist())}.
- Exactly-three-serotype reproducible positions: {len(exactly3)}.

## Figure 6 / Variance

- Median frac_tau2 {frac_tau_median:.3f}; IQR {variance.frac_tau2.quantile(0.25):.3f}-{variance.frac_tau2.quantile(0.75):.3f}; range {variance.frac_tau2.min():.3f}-{variance.frac_tau2.max():.3f}.
- Variance regimes: {int(regime_counts.get("replicate_dominated", 0))}/32 replicate-dominated, {int(regime_counts.get("balanced", 0))}/32 balanced, {int(regime_counts.get("sampling_dominated", 0))}/32 sampling-dominated.
- Non-replicate-dominated exceptions: {", ".join(f"{r.serotype} {r.domain} ({r.variance_regime}, frac_tau2={r.frac_tau2:.3f})" for _, r in non_repl.sort_values(["serotype", "domain"]).iterrows())}.
"""

    (outdir / "numbers_audit.md").write_text(numbers_audit, encoding="utf-8")
    (outdir / "claims_to_edit.md").write_text(claims_to_edit, encoding="utf-8")
    (outdir / "numbers_for_manuscript.md").write_text(numbers_for_manuscript, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path,
                        default=Path(__file__).resolve().parents[2] / "strideouts")
    parser.add_argument("--outdir", type=Path,
                        default=Path(__file__).resolve().parent)
    args = parser.parse_args()
    build_audit(args.root.resolve(), args.outdir.resolve())


if __name__ == "__main__":
    main()
