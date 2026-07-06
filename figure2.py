"""Figure 2 — Reproducibility emerges at the domain scale.

Biological message : per-residue effects are largely irreproducible, but coarsening
    to secondary-structure → motif → domain recovers reproducibility, which is
    therefore a domain-scale property; almost all licensing occurs at the domain scale.
Why this figure exists : it establishes the spatial scale at which the analysis is
    interpretable, motivating the domain-scale focus of the rest of the paper.
Generated from : outputs_s7/F7_rho_vs_scale_catalytic.parquet (ρ across scales) and
    outputs_s7/F2_resolution_census.parquet (gated-scale census).
Placement : MAIN. Merges the former scale-trajectory and resolution-census figures;
    the trivial chain/protein/complex scales (ρ→1 by construction) are dropped.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import (
    DOUBLE_COL,
    SCALE_ORDER,
    SCALE_SHORT,
    SEROTYPE_ORDER,
    Paths,
)
from utils import load_s7_figure, panel_letter, save_figure

# residue → secondary_structure → motif → domain (coarser scales are omitted:
# chain/protein/complex ρ saturate at ~1.0 for essentially every locus and carry
# no discriminating signal).
_KEEP_SCALES = ("residue", "secondary_structure", "motif", "domain")


def _panel_scale(ax, df: pd.DataFrame, serotypes: list[str]) -> None:
    keep_idx = [SCALE_ORDER.index(s) for s in _KEEP_SCALES]
    df = df[df["scale_index"].isin(keep_idx)]
    for sero in serotypes:
        sub = df[df["serotype"] == sero]
        color = styles.serotype_color(sero)
        for _, g in sub.groupby("canon_label"):
            g = g.sort_values("scale_index")
            ax.plot(g["scale_index"], g["rho"], color=color, lw=0.6, alpha=0.28,
                    zorder=2)
        mean_traj = sub.groupby("scale_index")["rho"].mean()
        ax.plot(mean_traj.index, mean_traj.values, color=color, lw=2.0, zorder=3,
                label=sero)
    ax.set_xticks(keep_idx)
    ax.set_xticklabels([SCALE_SHORT[SCALE_ORDER[i]] for i in keep_idx])
    ax.set_ylim(0, 1.02)
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_ylabel("reproducibility ρ")
    ax.set_xlabel("spatial scale (fine → coarse)")
    ax.set_title("ρ saturates by the domain scale", fontsize=styles.FS_LABEL, loc="right")


def _panel_census(ax, df: pd.DataFrame, serotypes: list[str]) -> None:
    df = df.copy()
    df["gated_scale_index"] = pd.to_numeric(df["gated_scale_index"], errors="coerce")
    df["n_loci"] = pd.to_numeric(df["n_loci"], errors="coerce").fillna(0)
    levels = (df[["gated_scale_index", "gated_scale_level"]].drop_duplicates()
              .sort_values("gated_scale_index"))
    level_names = list(levels["gated_scale_level"])
    ramp = plt.get_cmap(styles.SEQ_CMAP, max(2, len(level_names)))
    x = np.arange(len(serotypes))
    bottoms = np.zeros(len(serotypes))
    for k, lv in enumerate(level_names):
        h = np.array([df[(df["serotype"] == s) & (df["gated_scale_level"] == lv)]
                      ["n_loci"].sum() for s in serotypes], dtype=float)
        ax.bar(x, h, bottom=bottoms, width=0.66, color=ramp(k), edgecolor="white",
               linewidth=0.5, label=lv.replace("_", " "))
        bottoms += h
    ax.set_xticks(x)
    ax.set_xticklabels(serotypes, rotation=0)
    ax.set_ylabel("gated loci (count)")
    ax.set_title("loci gate at the domain scale", fontsize=styles.FS_LABEL, loc="right")
    ax.margins(y=0.08)
    ax.legend(title="gated scale", loc="upper center", bbox_to_anchor=(0.5, -0.16),
              ncol=min(3, len(level_names)), fontsize=styles.FS_ANNOT)


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    scale = load_s7_figure(paths, "F7_rho_vs_scale_catalytic").copy()
    scale["scale_index"] = pd.to_numeric(scale["scale_index"], errors="coerce")
    scale["rho"] = pd.to_numeric(scale["rho"], errors="coerce")
    census = load_s7_figure(paths, "F2_resolution_census")

    serotypes = [s for s in SEROTYPE_ORDER if s in set(scale["serotype"])]
    serotypes += sorted(set(scale["serotype"].astype(str)) - set(serotypes))

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 2.7),
        gridspec_kw={"width_ratios": [1.5, 1.0], "wspace": 0.28})
    _panel_scale(axa, scale, serotypes)
    axa.legend(title="serotype", loc="lower right", fontsize=styles.FS_ANNOT,
               ncol=1)
    panel_letter(axa, "A")
    _panel_census(axb, census, serotypes)
    panel_letter(axb, "B")

    fig.suptitle("Reproducibility emerges at the domain scale", x=0.5, y=1.02)
    return save_figure(fig, paths, "figure2_spatial_resolution")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
