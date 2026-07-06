"""Supplementary Figure S5 — Cross-serotype reproducibility matrix.

A compact 4×4 overview of which serotypes behave most similarly at the dynamic
level: each cell is the Pearson correlation between two serotypes' residue-level
reproducibility (ρ) profiles across shared canonical positions. Rows/columns are
ordered by hierarchical clustering and a dendrogram is drawn above. As the prompt
notes, this introduces no new statistics beyond pairwise correlation and clustering
*for visualization*.

Input: ``outputs_s7/F1_reproducibility_landscape.parquet`` (per-serotype, per-residue
ρ), pivoted to positions × serotypes.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import ONEHALF_COL, SEROTYPE_ORDER, Paths
from utils import load_s7_figure, save_figure


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F1_reproducibility_landscape").copy()
    df["rho_residue"] = pd.to_numeric(df["rho_residue"], errors="coerce")

    mat = df.pivot_table(index="canon_label", columns="serotype",
                         values="rho_residue", aggfunc="mean")
    serotypes = [s for s in SEROTYPE_ORDER if s in mat.columns]
    serotypes += [s for s in mat.columns if s not in serotypes]
    mat = mat[serotypes]

    # pairwise Pearson correlation of ρ profiles (pairwise-complete positions)
    corr = mat.corr(method="pearson", min_periods=2)

    fig = plt.figure(figsize=(ONEHALF_COL, ONEHALF_COL * 0.95))
    order = list(corr.columns)
    dendro_ok = False
    corr_vals = corr.to_numpy()
    finite = np.isfinite(corr_vals).all()
    if len(order) >= 3 and finite:
        dist = 1.0 - corr_vals
        np.fill_diagonal(dist, 0.0)
        dist = (dist + dist.T) / 2.0  # enforce symmetry for squareform
        degenerate = float(np.nanmax(dist)) < 1e-9  # all profiles identical
        if not degenerate:
            try:
                from scipy.cluster.hierarchy import dendrogram, linkage
                from scipy.spatial.distance import squareform

                z = linkage(squareform(dist, checks=False), method="average")
                gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 4.0], hspace=0.05)
                axd = fig.add_subplot(gs[0])
                dendro = dendrogram(
                    z, labels=list(corr.columns), no_labels=True,
                    color_threshold=0,
                    above_threshold_color=styles.OKABE_ITO["grey"], ax=axd)
                axd.set_axis_off()
                order = [corr.columns[i] for i in dendro["leaves"]]
                ax = fig.add_subplot(gs[1])
                dendro_ok = True
            except Exception:
                ax = fig.add_subplot(1, 1, 1)
        else:
            ax = fig.add_subplot(1, 1, 1)
    else:
        ax = fig.add_subplot(1, 1, 1)

    cm = corr.loc[order, order]
    im = ax.imshow(cm.to_numpy(dtype=float), cmap="cividis", vmin=-1, vmax=1,
                   aspect="equal")
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha="right")
    ax.set_yticklabels(order)
    for i in range(len(order)):
        for j in range(len(order)):
            v = cm.to_numpy()[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=styles.FS_ANNOT,
                        color="white" if v < 0.35 else "black")
    ax.set_title("Cross-serotype similarity of ρ profiles" if not dendro_ok
                 else None)
    if dendro_ok:
        fig.suptitle("Cross-serotype similarity of ρ profiles", y=0.98)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    cbar.set_label("Pearson r (ρ profiles)", fontsize=styles.FS_TICK)
    cbar.outline.set_linewidth(0.4)
    return save_figure(fig, paths, "Supplementary_Figure_S5")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
