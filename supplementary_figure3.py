"""Supplementary Figure S3 — Domain summary ranking.

A one-page summary of domain behaviour: for every (chain, domain) region the median
and mean ρ, mean coherence, and mean variance contribution (τ² fraction) are
aggregated across serotypes and shown as horizontal ranked bars (median ρ), with a
companion dot panel for coherence and variance contribution. Domains are ranked by
median ρ and coloured by reproducibility; catalytic domains are flagged.

Inputs (existing outputs; aggregation only, no recomputation):
``outputs_s7/T2_domain_rho_signed_effect.parquet`` (rho_domain, coherence_domain) and
``outputs_s7/T5_variance_component_budget.parquet`` (frac_tau2).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

import styles
from figure_config import CATALYTIC_DOMAINS, DOUBLE_COL, Paths
from utils import load_s7_figure, panel_letter, save_figure


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    t2 = load_s7_figure(paths, "T2_domain_rho_signed_effect").copy()
    t5 = load_s7_figure(paths, "T5_variance_component_budget").copy()
    for d, cols in ((t2, ["rho_domain", "coherence_domain"]), (t5, ["frac_tau2"])):
        d["chain"] = d["chain"].astype(str)
        d["domain"] = d["domain"].astype(str)
        for c in cols:
            d[c] = pd.to_numeric(d[c], errors="coerce")

    # aggregate across serotypes per (chain, domain)
    agg = (
        t2.groupby(["chain", "domain"])
        .agg(median_rho=("rho_domain", "median"),
             mean_rho=("rho_domain", "mean"),
             coherence=("coherence_domain", "mean"))
        .reset_index()
    )
    var = (
        t5.groupby(["chain", "domain"])["frac_tau2"].mean()
        .reset_index().rename(columns={"frac_tau2": "variance_contribution"})
    )
    agg = agg.merge(var, on=["chain", "domain"], how="left")
    agg = agg.sort_values("median_rho", ascending=True).reset_index(drop=True)

    labels = [
        (f"* {dom} ({chain})" if dom in CATALYTIC_DOMAINS else f"{dom} ({chain})")
        for chain, dom in zip(agg["chain"], agg["domain"])
    ]
    y = np.arange(len(agg))
    norm = Normalize(0, 1)
    cmap = plt.get_cmap(styles.SEQ_CMAP)

    fig, (ax, axc) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 0.34 * max(1, len(agg)) + 1.0),
        gridspec_kw={"width_ratios": [2.6, 1.4], "wspace": 0.06}, sharey=True,
    )

    # -- ranked median-ρ bars --------------------------------------------
    ax.barh(y, agg["median_rho"], color=cmap(norm(agg["median_rho"])),
            edgecolor="white", linewidth=0.5, zorder=2)
    # mean ρ as a thin marker on each bar
    ax.scatter(agg["mean_rho"], y, marker="|", s=90, color=styles.OKABE_ITO["black"],
               linewidths=1.0, zorder=3, label="mean ρ")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=styles.FS_ANNOT)
    for tick, dom in zip(ax.get_yticklabels(), agg["domain"], strict=True):
        if dom in CATALYTIC_DOMAINS:
            tick.set_color(styles.CATALYTIC_ACCENT)
            tick.set_fontweight("bold")
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("ρ (domain, median across serotypes)")
    ax.set_title("Domain reproducibility ranking")
    ax.legend(loc="lower right", fontsize=styles.FS_ANNOT)
    panel_letter(ax, "A")

    # -- companion: coherence + variance contribution --------------------
    axc.scatter(agg["coherence"], y, s=30, color=styles.COHERENT_COLOR,
                edgecolors="white", linewidths=0.4, label="coherence", zorder=2)
    axc.scatter(agg["variance_contribution"], y, s=30, marker="D",
                color=styles.VARIANCE_COLORS["tau2"], edgecolors="white",
                linewidths=0.4, label="τ² fraction", zorder=2)
    axc.set_xlim(0, 1.0)
    axc.set_xlabel("fraction")
    axc.set_title("coherence · variance")
    axc.legend(loc="lower right", fontsize=styles.FS_ANNOT)
    panel_letter(axc, "B")

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=[ax, axc], fraction=0.02, pad=0.02)
    cbar.set_label("median ρ", fontsize=styles.FS_TICK)
    cbar.outline.set_linewidth(0.4)
    return save_figure(fig, paths, "Supplementary_Figure_S3")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
