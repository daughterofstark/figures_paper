"""Figure 1 — Reproducibility landscape (redesign of S7 F1).

The original F1 is a spaghetti plot overlaying every serotype. Here it is redrawn
as **serotype-faceted small multiples**: one panel per serotype showing the
per-residue reproducibility (rho) along the sequence, with catalytic residues
highlighted and a rolling mean to reveal the trend. This removes overplotting and
makes per-serotype structure legible.

Consumes: ``outputs_s7/F1_reproducibility_landscape.parquet``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import (
    DOUBLE_COL,
    PROVISIONAL_RHO_STAR,
    SEROTYPE_ORDER,
    Paths,
)
from utils import (
    is_catalytic,
    load_s7_figure,
    panel_letter,
    residue_order,
    save_figure,
)

_LETTERS = "ABCDEFGH"


def _rolling(y: np.ndarray, window: int = 9) -> np.ndarray:
    if len(y) < window:
        return y
    kernel = np.ones(window) / window
    return np.convolve(y, kernel, mode="same")


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F1_reproducibility_landscape")

    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))
    n = max(1, len(serotypes))

    fig, axes = plt.subplots(
        n, 1, figsize=(DOUBLE_COL, 1.05 * n + 0.6), sharex=True, squeeze=False
    )
    axes = axes[:, 0]

    order = residue_order(df["canon_label"])
    xpos = {lab: i for i, lab in enumerate(order)}

    for i, sero in enumerate(serotypes):
        ax = axes[i]
        sub = df[df["serotype"] == sero].copy()
        sub["x"] = sub["canon_label"].map(xpos)
        sub = sub.sort_values("x")
        x = sub["x"].to_numpy()
        y = pd.to_numeric(sub["rho_residue"], errors="coerce").to_numpy()

        ax.axhline(
            PROVISIONAL_RHO_STAR, color=styles.OKABE_ITO["grey"],
            lw=0.6, ls=(0, (4, 3)), zorder=1,
        )
        ax.plot(
            x, y, color=styles.serotype_color(sero), lw=0.8, alpha=0.55,
            zorder=2, marker="o", ms=2.5, mfc=styles.serotype_color(sero),
            mec="none",
        )
        if len(x) >= 9:
            ax.plot(x, _rolling(y), color=styles.serotype_color(sero), lw=1.6,
                    zorder=3)
        # highlight catalytic residues
        cat = sub[is_catalytic(sub["domain"])]
        if not cat.empty:
            ax.scatter(
                cat["x"], pd.to_numeric(cat["rho_residue"], errors="coerce"),
                s=26, facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                linewidths=1.1, zorder=4,
            )
        ax.set_ylim(0, 1.02)
        ax.set_ylabel(f"{sero}\nρ", rotation=0, ha="right", va="center",
                      fontsize=styles.FS_TICK)
        ax.margins(x=0.01)
        if i == 0:
            panel_letter(ax, "A")

    axes[-1].set_xlabel("canonical residue position (ordered by chain, number)")
    # sparse x ticks with labels when few positions
    if len(order) <= 24:
        axes[-1].set_xticks(range(len(order)))
        axes[-1].set_xticklabels(order, rotation=90, fontsize=styles.FS_ANNOT)

    # shared legend for the catalytic marker + provisional line
    handles = [
        plt.Line2D([], [], marker="o", ls="none", mfc="none",
                   mec=styles.CATALYTIC_ACCENT, ms=6, mew=1.1,
                   label="catalytic residue"),
        plt.Line2D([], [], color=styles.OKABE_ITO["grey"], lw=0.8,
                   ls=(0, (4, 3)), label=f"provisional ρ* = {PROVISIONAL_RHO_STAR}"),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 1.0),
               ncol=2)
    fig.suptitle("Per-residue reproducibility landscape by serotype",
                 x=0.5, y=1.005)
    fig.tight_layout()
    return save_figure(fig, paths, "figure1_reproducibility_landscape")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
