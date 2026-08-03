"""Figure 1: observed item-level magnitude stability and coordinate presence.

Inputs are frozen S7 item-level rho rows and the frozen S5 presence/landmark
table. Absence remains an explicit state and is never converted to rho=0.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import ListedColormap, Normalize

import styles
from figure_config import DOUBLE_COL, SEROTYPE_ORDER, Paths
from phase3j_support import sha256, write_json
from utils import load_s7_figure, load_table, panel_letter, residue_order, save_figure

TITLE = "Observed item-level magnitude stability varies across coordinates"
PANEL_TITLES = (
    "Dataset-specific item rho and explicit absence",
    "Across-present median and presence count",
)
FIGURE_STRINGS = (TITLE, *PANEL_TITLES, "Item-level rho", "Across-present median rho")


def prepare_data(paths: Paths) -> tuple[pd.DataFrame, pd.DataFrame, list[str]]:
    """Validate and align frozen item rows and presence metadata."""
    items = load_s7_figure(paths, "F1_reproducibility_landscape").copy()
    presence = load_table(paths, "s5", "position_conservation").copy()
    required_items = {"serotype", "canon_label", "chain", "rho_residue"}
    required_presence = {
        "canon_label",
        "chain",
        "n_serotypes_present",
        "rho_residue_median",
        "is_catalytic_triad",
    }
    if missing := required_items - set(items.columns):
        raise ValueError(f"Figure 1 item input missing columns: {sorted(missing)}")
    if missing := required_presence - set(presence.columns):
        raise ValueError(f"Figure 1 presence input missing columns: {sorted(missing)}")
    if items.duplicated(["serotype", "canon_label"]).any():
        raise ValueError("Figure 1 item input has duplicate dataset-coordinate rows")
    if presence["canon_label"].duplicated().any():
        raise ValueError("Figure 1 presence input has duplicate coordinates")

    items["rho_residue"] = pd.to_numeric(items["rho_residue"], errors="raise")
    presence["n_serotypes_present"] = pd.to_numeric(
        presence["n_serotypes_present"], errors="raise"
    ).astype(int)
    order = residue_order(presence["canon_label"])
    if set(items["canon_label"].astype(str)) != set(order):
        raise ValueError("Figure 1 item and presence coordinate sets differ")
    present_counts = items.groupby("canon_label")["serotype"].nunique()
    expected_counts = presence.set_index("canon_label")["n_serotypes_present"]
    if not present_counts.reindex(order).equals(expected_counts.reindex(order)):
        raise ValueError("Figure 1 presence counts disagree with frozen item rows")
    medians = items.groupby("canon_label")["rho_residue"].median().reindex(order)
    stored = pd.to_numeric(
        presence.set_index("canon_label")["rho_residue_median"], errors="raise"
    ).reindex(order)
    # S5 stores the same across-present median rounded to six decimals.
    authoritative_shape = len(items) == 924 and len(order) == 248
    if authoritative_shape and not np.allclose(medians, stored, rtol=0, atol=5.1e-7):
        raise ValueError("Figure 1 across-present medians disagree with S5")
    presence = presence.copy()
    presence["_rho_across_present"] = presence["canon_label"].map(medians)
    return items, presence, order


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    items, presence, order = prepare_data(paths)
    serotypes = [s for s in SEROTYPE_ORDER if s in set(items["serotype"])]
    serotypes += sorted(set(items["serotype"].astype(str)) - set(serotypes))
    xmap = {label: i for i, label in enumerate(order)}
    matrix = (
        items.pivot(index="serotype", columns="canon_label", values="rho_residue")
        .reindex(index=serotypes, columns=order)
        .to_numpy(float)
    )

    fig = plt.figure(figsize=(DOUBLE_COL, 3.15))
    gs = fig.add_gridspec(2, 2, width_ratios=[28, 1], height_ratios=[1.05, 1.55],
                          hspace=0.38, wspace=0.05)
    axa = fig.add_subplot(gs[0, 0])
    axb = fig.add_subplot(gs[1, 0], sharex=axa)
    cax = fig.add_subplot(gs[:, 1])
    norm = Normalize(0, 1)
    cmap = plt.get_cmap(styles.SEQ_CMAP).with_extremes(bad="#e6e6e6")
    im = axa.imshow(matrix, aspect="auto", cmap=cmap, norm=norm,
                    interpolation="nearest")
    axa.set_yticks(range(len(serotypes)))
    axa.set_yticklabels(serotypes)
    axa.tick_params(axis="x", labelbottom=False, bottom=False)
    axa.set_ylabel("Frozen dataset")
    axa.set_title(PANEL_TITLES[0], loc="left", fontsize=styles.FS_LABEL)
    panel_letter(axa, "A", dx=-0.045)

    summary = presence.set_index("canon_label").reindex(order)
    med = pd.to_numeric(summary["_rho_across_present"], errors="raise").to_numpy()
    n_present = summary["n_serotypes_present"].to_numpy(int)
    x = np.arange(len(order))
    sizes = 7 + 4 * n_present
    axb.scatter(x, med, c=med, cmap=cmap, norm=norm, s=sizes,
                edgecolors="white", linewidths=0.25, zorder=3)
    cat_mask = summary["is_catalytic_triad"].astype(bool).to_numpy()
    if cat_mask.any():
        axb.scatter(x[cat_mask], med[cat_mask], marker="v", s=32,
                    facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                    linewidths=0.9, zorder=4)
    partial = n_present < len(serotypes)
    axb.scatter(x[partial], np.full(partial.sum(), -0.075), marker="s",
                s=8 + 3 * n_present[partial], c=n_present[partial],
                cmap=ListedColormap(["#d9d9d9", "#969696", "#525252"]),
                vmin=1, vmax=max(3, len(serotypes) - 1), edgecolors="none",
                clip_on=False, zorder=3)
    axb.set_ylim(-0.13, 1.03)
    axb.set_yticks([0, 0.5, 1.0])
    axb.set_ylabel("Across-present\nmedian rho")
    axb.set_xlabel("Canonical coordinate (ordered by chain, number)")
    axb.set_title(PANEL_TITLES[1], loc="left", fontsize=styles.FS_LABEL)
    panel_letter(axb, "B", dx=-0.045)

    chains = summary["chain"].astype(str).to_numpy()
    changes = np.flatnonzero(chains[1:] != chains[:-1]) + 1
    for boundary in changes:
        for ax in (axa, axb):
            ax.axvline(boundary - 0.5, color="white" if ax is axa else "#bdbdbd",
                       lw=0.8, zorder=5)
    starts = np.r_[0, changes]
    ends = np.r_[changes, len(order)]
    for start, end in zip(starts, ends, strict=True):
        axa.text((start + end - 1) / 2, -0.55, chains[start],
                 ha="center", va="bottom", fontsize=styles.FS_ANNOT,
                 fontweight="bold", clip_on=False)

    axa.text(0.995, 1.08, "Grey = absent (not rho=0)", transform=axa.transAxes,
             ha="right", va="bottom", fontsize=styles.FS_ANNOT, color="#6f6f6f")
    axb.text(
        0.995,
        1.05,
        "Open triangles: catalytic-triad landmarks; lower strip: presence count <4",
        transform=axb.transAxes,
        ha="right",
        va="bottom",
        fontsize=styles.FS_ANNOT,
        color="#6f6f6f",
    )
    cbar = fig.colorbar(im, cax=cax)
    cbar.set_label("Item-level rho", fontsize=styles.FS_TICK)
    cbar.set_ticks([0, 0.5, 1.0])
    cbar.outline.set_linewidth(0.4)
    fig.suptitle(TITLE, x=0.5, y=0.995)
    written = save_figure(fig, paths, "figure1_reproducibility_landscape")

    counts = items.groupby("serotype").size().reindex(serotypes)
    stats = (
        items.groupby("serotype")["rho_residue"]
        .agg(["min", "median", "max"])
        .reindex(serotypes)
    )
    audit = {
        "dataset_counts": {k: int(v) for k, v in counts.items()},
        "union_coordinate_count": len(order),
        "partial_presence_coordinates": int((n_present < len(serotypes)).sum()),
        "absent_dataset_coordinate_cells": int(np.isnan(matrix).sum()),
        "rho_by_dataset": {
            k: {name: float(value) for name, value in row.items()}
            for k, row in stats.iterrows()
        },
        "sources": {
            str(paths.s7("F1_reproducibility_landscape.parquet")): sha256(
                paths.s7("F1_reproducibility_landscape.parquet")
            ),
            str(paths.table("s5", "position_conservation.parquet")): sha256(
                paths.table("s5", "position_conservation.parquet")
            ),
        },
    }
    write_json(paths.output / "figure1_numeric_audit.json", audit)
    (paths.output / "figure1_caption_draft.md").write_text(
        "**Figure 1. Observed item-level magnitude stability varies across "
        "coordinates.** (A) Frozen item-level rho values are shown separately "
        "for each dataset; grey cells denote coordinates absent from that "
        "dataset and are not rho=0. (B) The median uses only datasets in which "
        "each coordinate is present, with presence count retained below the "
        "axis. Open triangles mark the three catalytic-triad coordinates as "
        "orientation landmarks. Rho is a support-specific magnitude-stability "
        "statistic and does not encode sign agreement or probability of truth.\n"
    )
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
