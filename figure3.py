"""Figure 3 — Domain × serotype reproducibility heatmap (redesign of S7 F3).

Keeps the heatmap but improves it: domains ordered by mean rho (most reproducible
at top), catalytic domains flagged in the row labels, per-cell rho annotations, a
perceptually-uniform colour map with a slim colour bar, and hairline cell borders.

Consumes: ``outputs_s7/F3_domain_serotype_rho_heatmap.parquet``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


import styles
from figure_config import CATALYTIC_DOMAINS, ONEHALF_COL, SEROTYPE_ORDER, Paths
from utils import load_s7_figure, save_figure


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F3_domain_serotype_rho_heatmap").copy()
    df["rho_domain"] = pd.to_numeric(df["rho_domain"], errors="coerce")

    df["chain"] = df["chain"].astype(str)
    df["domain"] = df["domain"].astype(str)
    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))

    # Row identity is (chain, domain): the S5 matrix is keyed on
    # (serotype, chain, domain), so a domain label can recur across chains and
    # must not be collapsed. (serotype, chain, domain) is unique, so the pivot
    # never actually aggregates — each cell is a single ρ.
    mat = df.pivot_table(index=["chain", "domain"], columns="serotype",
                         values="rho_domain", aggfunc="mean")
    mat = mat.reindex(columns=serotypes)
    # order rows by mean rho (descending → most reproducible on top)
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
    regions = list(mat.index)  # list of (chain, domain) tuples

    n_rows, n_cols = mat.shape
    fig, ax = plt.subplots(
        figsize=(ONEHALF_COL, 0.32 * max(1, n_rows) + 1.1)
    )
    cmap = plt.get_cmap(styles.SEQ_CMAP).with_extremes(bad="#f0f0f0")
    im = ax.imshow(mat.to_numpy(dtype=float), aspect="auto", cmap=cmap,
                   vmin=0, vmax=1)

    ax.set_xticks(range(n_cols))
    ax.set_xticklabels(mat.columns, rotation=0)
    ax.set_yticks(range(n_rows))
    ylabels = [
        (f"* {dom} ({chain})" if dom in CATALYTIC_DOMAINS else f"{dom} ({chain})")
        for (chain, dom) in regions
    ]
    ax.set_yticklabels(ylabels)
    for tick, (_chain, dom) in zip(ax.get_yticklabels(), regions, strict=True):
        if dom in CATALYTIC_DOMAINS:
            tick.set_color(styles.CATALYTIC_ACCENT)
            tick.set_fontweight("bold")

    # hairline grid + per-cell annotations
    ax.set_xticks(np.arange(-0.5, n_cols, 1), minor=True)
    ax.set_yticks(np.arange(-0.5, n_rows, 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.8)
    ax.tick_params(which="minor", length=0)
    values = mat.to_numpy(dtype=float)
    for i in range(n_rows):
        for j in range(n_cols):
            v = values[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=styles.FS_ANNOT,
                        color="white" if v < 0.62 else "black")

    ax.set_title("Domain-scale reproducibility (ρ) across serotypes")
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("ρ (domain)", fontsize=styles.FS_TICK)
    cbar.outline.set_linewidth(0.4)
    ax.text(0.0, 1.02, "* catalytic machinery", transform=ax.transAxes,
            fontsize=styles.FS_ANNOT, color=styles.CATALYTIC_ACCENT)
    fig.tight_layout()
    return save_figure(fig, paths, "figure3_domain_serotype_heatmap")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
