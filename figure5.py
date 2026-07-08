"""Figure 5 — Cross-serotype reproducibility is heterogeneous.

Biological message : reproducibility conservation varies across positions; catalytic
    triad residues occupy different conservation classes; whole-profile serotype
    correlations are weak; positions reproducible in more serotypes tend to have
    higher median residue-level rho.
Why this figure exists : it answers how reproducible features are conserved across
    DENV1-4 without treating positions as independent biological replicates.
Generated from : outputs_s5/position_conservation.parquet (conservation per position),
    outputs_s7/F1_reproducibility_landscape.parquet (per-serotype ρ profiles → 4×4
    similarity).
Placement : MAIN. Cross-serotype conservation audit; signed-effect detail belongs in
    Figure 4 and Supplementary Figure S2.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import (
    load_s7_figure,
    load_table,
    panel_letter,
    save_figure,
)


def _panel_conservation(ax, cons: pd.DataFrame) -> None:
    d = cons.copy()
    d["n_serotypes_reproducible"] = pd.to_numeric(d["n_serotypes_reproducible"],
                                                  errors="coerce")
    counts = (d.groupby("n_serotypes_reproducible", observed=True).size()
              .reindex(range(5), fill_value=0))
    colors = [
        styles.CONSERVATION_COLORS["reproducible_none"],
        styles.CONSERVATION_COLORS["reproducible_some"],
        styles.CONSERVATION_COLORS["reproducible_some"],
        styles.CONSERVATION_COLORS["reproducible_majority"],
        styles.CONSERVATION_COLORS["reproducible_all"],
    ]
    x = np.arange(5)
    ax.bar(x, counts.to_numpy(), color=colors, edgecolor="white", linewidth=0.5)
    for xi, c in zip(x, counts.to_numpy(), strict=True):
        ax.text(xi, c + max(counts.max() * 0.02, 1), str(int(c)), ha="center",
                va="bottom", fontsize=styles.FS_ANNOT)
    cat = d[d.get("is_catalytic_triad", pd.Series(False, index=d.index)).astype(bool)]
    offsets = {1: -0.13, 2: 0.0, 3: 0.13, 4: 0.0, 0: 0.0}
    for _, row in cat.iterrows():
        n = int(row["n_serotypes_reproducible"])
        y0 = max(2.0, counts.loc[n] * 0.10)
        ax.scatter(n + offsets.get(n, 0.0), y0, s=62, facecolors="none",
                   edgecolors=styles.CATALYTIC_ACCENT, linewidths=1.2, zorder=3)
        ax.text(n + offsets.get(n, 0.0), y0 + max(counts.max() * 0.035, 2),
                row["canon_label"], ha="center", va="bottom",
                fontsize=styles.FS_ANNOT, color=styles.CATALYTIC_ACCENT)
    ax.set_xticks(x)
    ax.set_xticklabels(["0", "1", "2", "3", "4"])
    ax.set_xlabel("serotypes with residue ρ ≥ provisional ρ*")
    ax.set_ylabel("positions (count)")
    ax.set_title("conservation classes span 0-4 serotypes", fontsize=styles.FS_LABEL)
    ax.margins(y=0.15)


def _panel_matrix(fig, cell, land: pd.DataFrame) -> None:
    land = land.copy()
    land["rho_residue"] = pd.to_numeric(land["rho_residue"], errors="coerce")
    mat = land.pivot_table(index="canon_label", columns="serotype",
                           values="rho_residue", aggfunc="mean")
    serotypes = [s for s in SEROTYPE_ORDER if s in mat.columns]
    serotypes += [s for s in mat.columns if s not in serotypes]
    mat = mat[serotypes]
    corr = mat.corr(method="pearson", min_periods=2)
    order = list(corr.columns)
    corr_vals = corr.to_numpy()

    axd = None
    if len(order) >= 3 and np.isfinite(corr_vals).all():
        dist = 1.0 - corr_vals
        np.fill_diagonal(dist, 0.0)
        dist = (dist + dist.T) / 2.0
        if float(np.nanmax(dist)) >= 1e-9:
            try:
                from scipy.cluster.hierarchy import dendrogram, linkage
                from scipy.spatial.distance import squareform

                z = linkage(squareform(dist, checks=False), method="average")
                sub = cell.subgridspec(2, 1, height_ratios=[1, 4], hspace=0.05)
                axd = fig.add_subplot(sub[0])
                dd = dendrogram(z, no_labels=True, color_threshold=0,
                                above_threshold_color="#9e9e9e", ax=axd)
                axd.set_axis_off()
                order = [corr.columns[i] for i in dd["leaves"]]
                ax = fig.add_subplot(sub[1])
            except Exception:
                ax = fig.add_subplot(cell)
        else:
            ax = fig.add_subplot(cell)
    else:
        ax = fig.add_subplot(cell)

    cm = corr.loc[order, order]
    im = ax.imshow(cm.to_numpy(float), cmap=styles.CORR_CMAP, vmin=-1, vmax=1,
                   aspect="equal")
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=styles.FS_ANNOT)
    ax.set_yticklabels(order, fontsize=styles.FS_ANNOT)
    for i in range(len(order)):
        for j in range(len(order)):
            v = cm.to_numpy()[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=styles.FS_ANNOT,
                        color="white" if v < 0.35 else "black")
    # title on the top-most axis (dendrogram if present) so nothing overlaps
    top_ax = axd if axd is not None else ax
    top_ax.set_title("serotype profile correlations are weak", fontsize=styles.FS_LABEL,
                     pad=4)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r", fontsize=styles.FS_TICK)
    cbar.set_ticks([-1, 0, 1])
    cbar.outline.set_linewidth(0.4)
    return axd if axd is not None else ax


def _panel_rho_by_conservation(ax, cons: pd.DataFrame) -> None:
    df = cons[["canon_label", "n_serotypes_reproducible", "rho_residue_median",
               "is_catalytic_triad"]].copy()
    df["n_serotypes_reproducible"] = pd.to_numeric(
        df["n_serotypes_reproducible"], errors="coerce")
    df["rho_residue_median"] = pd.to_numeric(df["rho_residue_median"], errors="coerce")
    df = df.dropna(subset=["n_serotypes_reproducible", "rho_residue_median"])
    df["n_serotypes_reproducible"] = df["n_serotypes_reproducible"].astype(int)
    if df.empty:
        ax.text(0.5, 0.5, "no conservation rows", ha="center", va="center",
                transform=ax.transAxes, color=styles.NONSIG_COLOR)
        return
    xs = np.arange(5)
    grouped = [df.loc[df["n_serotypes_reproducible"] == n, "rho_residue_median"]
               for n in xs]
    bp = ax.boxplot(grouped, positions=xs, widths=0.55, patch_artist=True,
                    showfliers=False, medianprops={"color": styles.OKABE_ITO["black"],
                                                   "linewidth": 1.0},
                    whiskerprops={"color": "#777777", "linewidth": 0.7},
                    capprops={"color": "#777777", "linewidth": 0.7})
    colors = [
        styles.CONSERVATION_COLORS["reproducible_none"],
        styles.CONSERVATION_COLORS["reproducible_some"],
        styles.CONSERVATION_COLORS["reproducible_some"],
        styles.CONSERVATION_COLORS["reproducible_majority"],
        styles.CONSERVATION_COLORS["reproducible_all"],
    ]
    for patch, color in zip(bp["boxes"], colors, strict=True):
        patch.set_facecolor(color)
        patch.set_alpha(0.62)
        patch.set_edgecolor("#777777")
        patch.set_linewidth(0.6)
    offsets = ((np.arange(len(df)) % 9) - 4) * 0.018
    xj = df["n_serotypes_reproducible"].to_numpy(float) + offsets
    cat = df["is_catalytic_triad"].astype(bool).to_numpy()
    ax.scatter(xj[~cat], df["rho_residue_median"].to_numpy()[~cat], s=9,
               c="#9e9e9e", alpha=0.34, edgecolors="none", zorder=2)
    ax.scatter(xj[cat], df["rho_residue_median"].to_numpy()[cat], s=52,
               facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
               linewidths=1.1, zorder=4, label="catalytic triad")
    for _, row in df[df["is_catalytic_triad"].astype(bool)].iterrows():
        ax.text(row["n_serotypes_reproducible"] + 0.08, row["rho_residue_median"],
                row["canon_label"], ha="left", va="center",
                fontsize=styles.FS_ANNOT, color=styles.CATALYTIC_ACCENT)
    ax.set_xlim(-0.55, 4.55)
    ax.set_ylim(0, 1.0)
    ax.set_xticks(xs)
    ax.set_xlabel("serotypes with residue ρ ≥ provisional ρ*")
    ax.set_ylabel("median residue-level ρ")
    ax.set_title("more conserved positions have higher median ρ", fontsize=styles.FS_LABEL)
    ax.legend(loc="upper left", fontsize=styles.FS_ANNOT)


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    cons = load_table(paths, "s5", "position_conservation")
    land = load_s7_figure(paths, "F1_reproducibility_landscape")

    fig = plt.figure(figsize=(DOUBLE_COL, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.35, 1.0],
                          hspace=0.42, wspace=0.34)
    axa = fig.add_subplot(gs[0, 0])
    _panel_conservation(axa, cons)
    panel_letter(axa, "A")

    top_b = _panel_matrix(fig, gs[0, 1], land)
    panel_letter(top_b, "B")

    axc = fig.add_subplot(gs[1, :])
    _panel_rho_by_conservation(axc, cons)
    panel_letter(axc, "C")

    fig.suptitle("Cross-serotype reproducibility is heterogeneous across positions",
                 x=0.5, y=1.0)
    return save_figure(fig, paths, "figure5_cross_serotype")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
