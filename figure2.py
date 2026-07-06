"""Figure 2 — Achieved-resolution census (redesign of S7 F2).

Keeps the stacked-bar concept (the design's intent) but improves typography,
spacing, ordering (coarsest scale at the bottom), a sequential colour ramp keyed
to spatial scale, and in-bar count labels.

Consumes: ``outputs_s7/F2_resolution_census.parquet``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors

import styles
from figure_config import SEROTYPE_ORDER, SINGLE_COL, Paths
from utils import load_s7_figure, save_figure


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F2_resolution_census").copy()
    df["gated_scale_index"] = pd.to_numeric(df["gated_scale_index"], errors="coerce")
    df["n_loci"] = pd.to_numeric(df["n_loci"], errors="coerce").fillna(0)

    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))

    # scale levels present, ordered fine→coarse (by gated_scale_index)
    levels = (
        df[["gated_scale_index", "gated_scale_level"]]
        .drop_duplicates()
        .sort_values("gated_scale_index")
    )
    level_names = list(levels["gated_scale_level"])
    ramp = plt.get_cmap(styles.SEQ_CMAP, max(2, len(level_names)))
    level_color = {lv: ramp(i) for i, lv in enumerate(level_names)}

    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.82))
    x = np.arange(len(serotypes))
    bottoms = np.zeros(len(serotypes))
    for lv in level_names:
        heights = np.array(
            [
                df[(df["serotype"] == s) & (df["gated_scale_level"] == lv)]["n_loci"].sum()
                for s in serotypes
            ],
            dtype=float,
        )
        ax.bar(x, heights, bottom=bottoms, width=0.66, color=level_color[lv],
               edgecolor="white", linewidth=0.5, label=lv.replace("_", " "))
        for xi, (h, b) in enumerate(zip(heights, bottoms, strict=True)):
            if h > 0:
                ax.text(xi, b + h / 2, f"{int(h)}", ha="center", va="center",
                        fontsize=styles.FS_ANNOT,
                        color=_readable_text(level_color[lv]))
        bottoms += heights

    ax.set_xticks(x)
    ax.set_xticklabels(serotypes)
    ax.set_ylabel("gated loci (count)")
    ax.set_xlabel("serotype")
    ax.set_title("Achieved-resolution census")
    ax.margins(y=0.08)
    ax.legend(title="gated spatial scale", loc="upper center",
              bbox_to_anchor=(0.5, -0.16), ncol=min(3, len(level_names)))
    fig.tight_layout()
    return save_figure(fig, paths, "figure2_resolution_census")


def _readable_text(rgba: tuple[float, float, float, float]) -> str:
    r, g, b = colors.to_rgb(rgba)
    luminance = 0.299 * r + 0.587 * g + 0.114 * b
    return "white" if luminance < 0.5 else "black"


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
