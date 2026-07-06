"""Figure 7 — ρ-vs-scale trajectories for catalytic regions (redesign of S7 F7).

Keeps the concept (reproducibility as a function of spatial scale) but reduces
clutter: serotypes become small-multiple panels, each catalytic residue is one
faint trajectory, and the panel's mean trajectory is drawn boldly on top. The
x-axis uses readable scale abbreviations (res → cplx).

Consumes: ``outputs_s7/F7_rho_vs_scale_catalytic.parquet``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import (
    DOUBLE_COL,
    PROVISIONAL_RHO_STAR,
    SCALE_ORDER,
    SCALE_SHORT,
    SEROTYPE_ORDER,
    Paths,
)
from utils import load_s7_figure, panel_letter, save_figure

_LETTERS = "ABCDEFGH"


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F7_rho_vs_scale_catalytic").copy()
    df["scale_index"] = pd.to_numeric(df["scale_index"], errors="coerce")
    df["rho"] = pd.to_numeric(df["rho"], errors="coerce")

    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))
    n = max(1, len(serotypes))
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))

    fig, axes = plt.subplots(
        nrow, ncol, figsize=(DOUBLE_COL, 1.05 * nrow + 0.7),
        squeeze=False, sharex=True, sharey=True,
    )
    flat = axes.ravel()
    xticks = sorted(df["scale_index"].dropna().unique())

    for k, sero in enumerate(serotypes):
        ax = flat[k]
        sub = df[df["serotype"] == sero]
        ax.axhline(PROVISIONAL_RHO_STAR, color=styles.OKABE_ITO["grey"], lw=0.6,
                   ls=(0, (4, 3)), zorder=1)
        for _, g in sub.groupby("canon_label"):
            g = g.sort_values("scale_index")
            ax.plot(g["scale_index"], g["rho"], color=styles.serotype_color(sero),
                    lw=0.7, alpha=0.4, zorder=2)
        # bold mean trajectory
        mean_traj = sub.groupby("scale_index")["rho"].mean()
        ax.plot(mean_traj.index, mean_traj.values,
                color=styles.serotype_color(sero), lw=2.0, zorder=3)
        ax.set_title(sero, fontsize=styles.FS_LABEL)
        ax.set_ylim(0, 1.02)
        panel_letter(ax, _LETTERS[k])

    for k in range(len(serotypes), len(flat)):
        flat[k].set_visible(False)

    for ax in flat[: len(serotypes)]:
        ax.set_xticks(xticks)
        ax.set_xticklabels(
            [SCALE_SHORT.get(SCALE_ORDER[int(i)], str(int(i)))
             if 0 <= int(i) < len(SCALE_ORDER) else str(int(i)) for i in xticks],
            rotation=0,
        )
    for k in range(0, len(serotypes), ncol):
        flat[k].set_ylabel("ρ")
    for ax in flat[max(0, len(serotypes) - ncol): len(serotypes)]:
        ax.set_xlabel("spatial scale (fine → coarse)")

    fig.suptitle("Reproducibility across spatial scale — catalytic regions", y=1.02)
    fig.tight_layout()
    return save_figure(fig, paths, "figure7_rho_vs_scale_catalytic")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
