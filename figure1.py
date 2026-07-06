"""Figure 1 — Reproducible dynamics localize to the catalytic core.

Biological message : along the NS2B–NS3 structure, per-position reproducibility (ρ)
    is highest at the catalytic machinery; reproducible dynamics are not uniform but
    localize to specific structural modules.
Why this figure exists : it is the manuscript's opening orientation figure — it maps
    *where* reproducible dynamics live, in structural context (chain + domain).
Generated from : outputs_s5/position_conservation.parquet (chain, domain,
    rho_residue_median, is_catalytic_triad; residue number parsed from canon_label).
Placement : MAIN. The structural, across-serotype "where" view. Per-serotype detail
    lives in Supplementary S1 (a genuinely different message).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize
from matplotlib.patches import Rectangle

import styles
from figure_config import CATALYTIC_DOMAINS, DOUBLE_COL, Paths
from utils import chain_bands, load_table, panel_letter, parse_canon_label, save_figure


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_table(paths, "s5", "position_conservation").copy()
    df["chain"] = df["chain"].astype(str)
    df["domain"] = df["domain"].astype(str)
    df["rho_residue_median"] = pd.to_numeric(df["rho_residue_median"], errors="coerce")
    df["_resid"] = df["canon_label"].map(lambda s: parse_canon_label(s)[1])
    df = df.sort_values(["chain", "_resid"]).reset_index(drop=True)
    df["_x"] = np.arange(len(df))

    fig = plt.figure(figsize=(DOUBLE_COL, 2.7))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    axd = fig.add_subplot(gs[1], sharex=ax)

    bands = chain_bands(list(df["chain"]))
    shades = ["#f4f4f4", "#eef2f7", "#f5f2ec", "#eef4ee"]
    band_shade = {c: shades[i % len(shades)] for i, c in enumerate(dict.fromkeys(df["chain"]))}
    for chain, s, e in bands:
        for a in (ax, axd):
            a.axvspan(s - 0.5, e - 0.5, color=band_shade[chain], zorder=0, lw=0)
        ax.text((s + e - 1) / 2, 1.04, chain, ha="center", va="bottom",
                fontsize=styles.FS_LABEL, fontweight="bold")

    norm = Normalize(0, 1)
    cmap = plt.get_cmap(styles.SEQ_CMAP)
    rho = df["rho_residue_median"].to_numpy()
    ax.vlines(df["_x"], 0, rho, color="#c8c8c8", lw=0.5, zorder=2)
    ax.scatter(df["_x"], rho, c=rho, cmap=cmap, norm=norm, s=16, zorder=3,
               edgecolors="white", linewidths=0.3)
    cat = df[df.get("is_catalytic_triad", pd.Series(False, index=df.index)).astype(bool)]
    if not cat.empty:
        ax.scatter(cat["_x"], cat["rho_residue_median"], s=52, marker="o",
                   facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                   linewidths=1.2, zorder=4, label="catalytic residue")
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0))
    ax.set_ylim(0, 1.10)
    ax.set_yticks([0.5, 1.0])
    ax.set_ylabel("reproducibility ρ\n(median over serotypes)")
    ax.tick_params(labelbottom=False)
    panel_letter(ax, "A")

    for (chain, dom), g in df.groupby(["chain", "domain"], sort=False):
        x0, x1 = g["_x"].min(), g["_x"].max()
        is_cat = dom in CATALYTIC_DOMAINS
        axd.add_patch(Rectangle(
            (x0 - 0.5, 0.15), (x1 - x0) + 1.0, 0.7,
            facecolor=(styles.CATALYTIC_ACCENT if is_cat else "#9ecae1"),
            alpha=0.9 if is_cat else 0.55,
            edgecolor=styles.OKABE_ITO["black"], linewidth=0.4))
        axd.text((x0 + x1) / 2, 0.5, dom, ha="center", va="center",
                 fontsize=styles.FS_ANNOT,
                 color="white" if is_cat else styles.OKABE_ITO["black"],
                 rotation=0 if (x1 - x0) >= 2 else 90)
    axd.set_ylim(0, 1)
    axd.set_yticks([])
    axd.set_ylabel("domains", rotation=0, ha="right", va="center",
                   fontsize=styles.FS_TICK)
    axd.set_xlabel("canonical residue position (ordered by chain, number)")
    panel_letter(axd, "B")
    if len(df) <= 30:
        axd.set_xticks(df["_x"])
        axd.set_xticklabels(df["canon_label"], rotation=90, fontsize=styles.FS_ANNOT)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=[ax, axd], fraction=0.02, pad=0.01)
    cbar.set_label("ρ", fontsize=styles.FS_TICK)
    cbar.set_ticks([0, 0.5, 1.0])
    cbar.outline.set_linewidth(0.4)
    fig.suptitle("Reproducible dynamics localize to the catalytic core of NS2B–NS3",
                 x=0.5, y=1.02)
    return save_figure(fig, paths, "figure1_reproducibility_landscape")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
