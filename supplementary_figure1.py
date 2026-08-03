"""Supplementary Figure S1: frozen per-dataset item-level rho profiles."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

import styles
from figure_config import DOUBLE_COL, PROVISIONAL_RHO_STAR, SEROTYPE_ORDER, Paths
from phase3j_support import sha256, write_json
from utils import load_s7_figure, load_table, panel_letter, residue_order, save_figure

TITLE = "Per-dataset item-level magnitude stability profiles"
REFERENCE_LABEL = f"provisional rho reference = {PROVISIONAL_RHO_STAR}"
FIGURE_STRINGS = (TITLE, REFERENCE_LABEL, "Item-level rho")


def prepare_data(paths: Paths) -> tuple[pd.DataFrame, list[str], set[str]]:
    """Validate frozen item rows and return deterministic coordinate order."""
    df = load_s7_figure(paths, "F1_reproducibility_landscape").copy()
    required = {"serotype", "canon_label", "chain", "rho_residue"}
    if missing := required - set(df.columns):
        raise ValueError(f"S1 input missing columns: {sorted(missing)}")
    if df.duplicated(["serotype", "canon_label"]).any():
        raise ValueError("S1 input has duplicate dataset-coordinate rows")
    df["rho_residue"] = pd.to_numeric(df["rho_residue"], errors="raise")
    order = residue_order(df["canon_label"])
    position = load_table(paths, "s5", "position_conservation")
    required_pos = {"canon_label", "is_catalytic_triad"}
    if missing := required_pos - set(position.columns):
        raise ValueError(f"S1 landmark input missing columns: {sorted(missing)}")
    catalytic = set(
        position.loc[position["is_catalytic_triad"].astype(bool), "canon_label"].astype(str)
    )
    return df, order, catalytic


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df, order, catalytic = prepare_data(paths)
    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))
    n = max(1, len(serotypes))

    fig, axes = plt.subplots(n, 1, figsize=(DOUBLE_COL, 1.0 * n + 0.6), sharex=True,
                             squeeze=False)
    axes = axes[:, 0]
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
        ax.scatter(x, y, color=color, s=9, alpha=0.72, edgecolors="white",
                   linewidths=0.2, zorder=2)
        cat = sub[sub["canon_label"].astype(str).isin(catalytic)]
        if not cat.empty:
            ax.scatter(cat["x"], cat["rho_residue"], s=28, marker="v",
                       facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                       linewidths=0.8, zorder=4)
        ax.set_ylim(0, 1.02)
        ax.set_yticks([0, 0.5, 1.0])
        ax.set_ylabel(f"{sero}\nrho", rotation=0, ha="right", va="center",
                      fontsize=styles.FS_TICK)
        ax.margins(x=0.01)
        if i == 0:
            panel_letter(ax, "A")

    axes[-1].set_xlabel("Canonical coordinate (ordered by chain, number)")
    if len(order) <= 30:
        axes[-1].set_xticks(range(len(order)))
        axes[-1].set_xticklabels(order, rotation=90, fontsize=styles.FS_ANNOT)

    handles = [
        plt.Line2D([], [], marker="v", ls="none", mfc="none",
                   mec=styles.CATALYTIC_ACCENT, ms=5, mew=0.8,
                   label="Catalytic-triad landmark"),
        plt.Line2D([], [], color="#d9d9d9", lw=0.8, ls=(0, (4, 3)),
                   label=REFERENCE_LABEL),
    ]
    fig.legend(handles=handles, loc="upper right", bbox_to_anchor=(0.995, 1.0), ncol=2,
               fontsize=styles.FS_ANNOT)
    fig.suptitle(TITLE, x=0.5, y=1.005)
    fig.tight_layout()
    written = save_figure(fig, paths, "Supplementary_Figure_S1")
    stats = df.groupby("serotype")["rho_residue"].agg(["count", "min", "median", "max"])
    write_json(
        paths.output / "Supplementary_Figure_S1_numeric_audit.json",
        {
            "no_interpolation": True,
            "no_smoothing": True,
            "provisional_reference": PROVISIONAL_RHO_STAR,
            "profiles": {
                k: {name: float(value) for name, value in row.items()}
                for k, row in stats.iterrows()
            },
            "source": str(paths.s7("F1_reproducibility_landscape.parquet")),
            "source_sha256": sha256(paths.s7("F1_reproducibility_landscape.parquet")),
        },
    )
    (paths.output / "Supplementary_Figure_S1_caption_draft.md").write_text(
        "**Supplementary Figure S1. Per-dataset item-level magnitude stability "
        "profiles.** Each point is one frozen item-level rho value. Coordinates "
        "absent from a dataset are left as gaps; no interpolation or smoothing "
        "is applied. The dashed line is the provisional rho reference of 0.5 "
        "and is not a calibrated decision threshold. Open triangles are "
        "catalytic-triad orientation landmarks.\n"
    )
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
