"""Supplementary Figure S1 — Structural domain map of NS2B–NS3.

A schematic of the canonical NS2B–NS3 sequence: residues are ordered along the
canonical axis, chains are shaded separately, every annotated domain is drawn as a
labelled block spanning the residues assigned to it, and per-position
reproducibility (ρ) is overlaid as a colour-coded lollipop. Catalytic-triad
residues are ringed. It answers, at a glance, *where along the protein the
reproducible regions sit*.

Every element is derived from existing canonical annotations
(``outputs_s5/position_conservation.parquet``): chain, domain, catalytic flag, the
median ρ per position, and the residue number parsed from ``canon_label``. Domain
spans are the min–max canonical residue of each (chain, domain) — a deterministic
aggregation of existing annotations, not manual coordinates.
"""
from __future__ import annotations

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

    fig = plt.figure(figsize=(DOUBLE_COL, 2.8))
    gs = fig.add_gridspec(2, 1, height_ratios=[3.0, 1.0], hspace=0.08)
    ax = fig.add_subplot(gs[0])
    axd = fig.add_subplot(gs[1], sharex=ax)

    # chain background bands
    bands = chain_bands(list(df["chain"]))
    band_shade = {c: shade for c, shade in zip(
        dict.fromkeys(df["chain"]), ["#f2f2f2", "#e6eef5", "#f5efe6", "#eef5ee"])}
    for chain, s, e in bands:
        for a in (ax, axd):
            a.axvspan(s - 0.5, e - 0.5, color=band_shade.get(chain, "#f2f2f2"),
                      zorder=0, lw=0)
        ax.text((s + e - 1) / 2, 1.06, chain, ha="center", va="bottom",
                fontsize=styles.FS_LABEL, fontweight="bold",
                color=styles.OKABE_ITO["black"])

    # ρ lollipops
    norm = Normalize(0, 1)
    cmap = plt.get_cmap(styles.SEQ_CMAP)
    rho = df["rho_residue_median"].to_numpy()
    ax.vlines(df["_x"], 0, rho, color=styles.OKABE_ITO["grey"], lw=0.5, zorder=2)
    ax.scatter(df["_x"], rho, c=rho, cmap=cmap, norm=norm, s=16, zorder=3,
               edgecolors="white", linewidths=0.3)
    cat = df[df.get("is_catalytic_triad", pd.Series(False, index=df.index)).astype(bool)]
    if not cat.empty:
        ax.scatter(cat["_x"], cat["rho_residue_median"], s=52, marker="o",
                   facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                   linewidths=1.2, zorder=4, label="catalytic triad")
        ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.0))
    ax.set_ylim(0, 1.12)
    ax.set_ylabel("ρ (median across serotypes)")
    ax.tick_params(labelbottom=False)
    ax.set_title("Structural domain map of NS2B–NS3 with per-position reproducibility")
    panel_letter(ax, "A")

    # domain ribbon
    for (chain, dom), g in df.groupby(["chain", "domain"], sort=False):
        x0, x1 = g["_x"].min(), g["_x"].max()
        is_cat = dom in CATALYTIC_DOMAINS
        axd.add_patch(Rectangle(
            (x0 - 0.5, 0.15), (x1 - x0) + 1.0, 0.7,
            facecolor=(styles.CATALYTIC_ACCENT if is_cat else styles.OKABE_ITO["sky"]),
            alpha=0.85 if is_cat else 0.55,
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

    # sparse x tick labels when few positions
    if len(df) <= 30:
        axd.set_xticks(df["_x"])
        axd.set_xticklabels(df["canon_label"], rotation=90, fontsize=styles.FS_ANNOT)

    sm = plt.cm.ScalarMappable(norm=norm, cmap=cmap)
    cbar = fig.colorbar(sm, ax=[ax, axd], fraction=0.02, pad=0.01)
    cbar.set_label("ρ", fontsize=styles.FS_TICK)
    cbar.outline.set_linewidth(0.4)
    return save_figure(fig, paths, "Supplementary_Figure_S1")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
