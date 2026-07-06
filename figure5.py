"""Figure 5 — Cross-serotype conservation (redesign of S7 F5).

The original conservation plot is unreadable. This redraws it as a **ranked dot
plot with a grouped-summary companion**: shared positions are ranked by the number
of serotypes in which they reproduce, dot size encodes the signed-and-reproducible
count, colour encodes conservation class, and catalytic residues are ringed. A
compact bar at the right summarises how many positions fall in each class.

Consumes: ``outputs_s5/position_conservation.parquet`` (richer than the S7 F5 slice).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import CONSERVATION_ORDER, DOUBLE_COL, Paths
from utils import load_table, panel_letter, parse_canon_label, save_figure


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_table(paths, "s5", "position_conservation").copy()
    for c in ("n_serotypes_reproducible", "n_serotypes_signed_reproducible",
              "frac_reproducible"):
        df[c] = pd.to_numeric(df[c], errors="coerce")

    # rank: most-reproducible first, then by shared fraction and residue order
    df["_key"] = df["canon_label"].map(parse_canon_label)
    df = df.sort_values(
        ["n_serotypes_reproducible", "frac_reproducible", "_key"],
        ascending=[True, True, True],
    ).reset_index(drop=True)

    fig, (ax, axb) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, max(2.2, 0.16 * len(df) + 1.0)),
        gridspec_kw={"width_ratios": [3.4, 1.0], "wspace": 0.32},
    )

    y = np.arange(len(df))
    class_color = [styles.CONSERVATION_COLORS.get(c, styles.OKABE_ITO["grey"])
                   for c in df["conservation_class"]]
    sizes = 18 + 26 * df["n_serotypes_signed_reproducible"].fillna(0)

    ax.hlines(y, 0, df["n_serotypes_reproducible"], color=styles.OKABE_ITO["grey"],
              lw=0.5, alpha=0.6, zorder=1)
    ax.scatter(df["n_serotypes_reproducible"], y, s=sizes, c=class_color,
               edgecolors="white", linewidths=0.4, zorder=2)
    cat = df[df.get("is_catalytic_triad", pd.Series(False, index=df.index)).astype(bool)]
    if not cat.empty:
        ax.scatter(cat["n_serotypes_reproducible"], cat.index, s=70,
                   facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                   linewidths=1.2, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(df["canon_label"], fontsize=styles.FS_ANNOT)
    ax.set_xlabel("serotypes in which the position reproduces (n)")
    ax.set_title("Cross-serotype conservation of reproducible positions")
    ax.margins(y=0.02)
    ax.set_xlim(-0.2, max(4, df["n_serotypes_reproducible"].max() or 4) + 0.4)
    panel_letter(ax, "A")

    # -- companion: count of positions per conservation class --------------
    present = [c for c in CONSERVATION_ORDER if c in set(df["conservation_class"])]
    counts = [int((df["conservation_class"] == c).sum()) for c in present]
    axb.barh(range(len(present)), counts,
             color=[styles.CONSERVATION_COLORS[c] for c in present],
             edgecolor="white", linewidth=0.5)
    axb.set_yticks(range(len(present)))
    axb.set_yticklabels([c.replace("reproducible_", "") for c in present],
                        fontsize=styles.FS_ANNOT)
    axb.invert_yaxis()
    for i, c in enumerate(counts):
        axb.text(c, i, f" {c}", va="center", ha="left", fontsize=styles.FS_ANNOT)
    axb.set_xlabel("positions (n)")
    axb.set_title("by class")
    panel_letter(axb, "B")

    # legends
    size_handles = [
        plt.Line2D([], [], marker="o", ls="none",
                   markerfacecolor=styles.OKABE_ITO["grey"],
                   markeredgecolor="white",
                   markersize=np.sqrt(18 + 26 * k),
                   label=f"{k} signed-reproducible")
        for k in (0, 2, 4)
    ]
    size_handles.append(
        plt.Line2D([], [], marker="o", ls="none", mfc="none",
                   mec=styles.CATALYTIC_ACCENT, ms=8, mew=1.2,
                   label="catalytic triad")
    )
    ax.legend(handles=size_handles, loc="lower right", bbox_to_anchor=(1.0, 0.0),
              fontsize=styles.FS_ANNOT)
    return save_figure(fig, paths, "figure5_cross_serotype_conservation")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
