"""Supplementary Figure S1 — Per-serotype residue-level reproducibility landscapes.

Biological message : the reproducibility landscape summarised in Main Figure 1 is
    shown separately for each serotype. Catalytic residues are functional landmarks,
    but the highest residue-level ρ values are distributed across the sequence.
Why this figure exists : it is the per-serotype detail behind Main Figure 1 (which
    shows the across-serotype median in structural context). Different message: here
    the emphasis is per-serotype reproducibility rather than structural localisation.
Generated from : outputs_s7/F1_reproducibility_landscape.parquet.
Placement : SUPPLEMENTARY (granular backup to Main Figure 1).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import DOUBLE_COL, PROVISIONAL_RHO_STAR, SEROTYPE_ORDER, Paths
from utils import is_catalytic, load_s7_figure, panel_letter, residue_order, save_figure


def _rolling(y: np.ndarray, window: int = 9) -> np.ndarray:
    if len(y) < window:
        return y
    return np.convolve(y, np.ones(window) / window, mode="same")


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F1_reproducibility_landscape").copy()
    df["rho_residue"] = pd.to_numeric(df["rho_residue"], errors="coerce")
    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))
    n = max(1, len(serotypes))

    fig, axes = plt.subplots(n, 1, figsize=(DOUBLE_COL, 1.0 * n + 0.6), sharex=True,
                             squeeze=False)
    axes = axes[:, 0]
    order = residue_order(df["canon_label"])
    xpos = {lab: i for i, lab in enumerate(order)}

    for i, sero in enumerate(serotypes):
        ax = axes[i]
        sub = df[df["serotype"] == sero].copy()
        sub["x"] = sub["canon_label"].map(xpos)
        sub = sub.sort_values("x")
        x = sub["x"].to_numpy()
        y = sub["rho_residue"].to_numpy()
        color = styles.serotype_color(sero)
        ax.axhline(PROVISIONAL_RHO_STAR, color="#d9d9d9", lw=0.6, ls=(0, (4, 3)),
                   zorder=1)
        ax.plot(x, y, color=color, lw=0.8, alpha=0.5, marker="o", ms=2.5,
                mfc=color, mec="none", zorder=2)
        if len(x) >= 9:
            ax.plot(x, _rolling(y), color=color, lw=1.6, zorder=3)
        cat = sub[is_catalytic(sub["domain"])]
        if not cat.empty:
            ax.scatter(cat["x"], cat["rho_residue"], s=26, facecolors="none",
                       edgecolors=styles.CATALYTIC_ACCENT, linewidths=1.1, zorder=4)
        ax.set_ylim(0, 1.02)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_ylabel(f"{sero}\nρ", rotation=0, ha="right", va="center",
                      fontsize=styles.FS_TICK)
        ax.margins(x=0.01)
        if i == 0:
            panel_letter(ax, "A")

    axes[-1].set_xlabel("canonical residue position (ordered by chain, number)")
    if len(order) <= 30:
        axes[-1].set_xticks(range(len(order)))
        axes[-1].set_xticklabels(order, rotation=90, fontsize=styles.FS_ANNOT)

    handles = [
        plt.Line2D([], [], marker="o", ls="none", mfc="none",
                   mec=styles.CATALYTIC_ACCENT, ms=6, mew=1.1, label="catalytic residue"),
        plt.Line2D([], [], color="#d9d9d9", lw=0.8, ls=(0, (4, 3)),
                   label=f"provisional ρ* = {PROVISIONAL_RHO_STAR}"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 1.0), ncol=2,
               fontsize=styles.FS_ANNOT)
    fig.suptitle("Per-serotype residue-level ρ landscapes are heterogeneous", x=0.5,
                 y=1.005)
    fig.tight_layout()
    return save_figure(fig, paths, "Supplementary_Figure_S1")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
