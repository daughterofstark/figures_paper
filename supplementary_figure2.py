"""Supplementary Figure S2 — Signed dynamic effects along the NS2B–NS3 sequence.

Biological message : the significant signed effects of Main Figure 4, placed on the
    sequence axis, are distributed across NS2B and NS3; catalytic residues are
    highlighted as functional reference points.
Why this figure exists : it is the positional (per-residue) companion to the Figure 4
    volcano — same effects, arranged along the sequence for structural readers.
Generated from : outputs_s4/significance_screen.parquet.
Placement : SUPPLEMENTARY (positional backup to Main Figure 4).
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

import styles
from figure_config import DOUBLE_COL, Paths
from utils import (
    chain_bands,
    is_catalytic,
    load_table,
    objective_labels,
    parse_canon_label,
    save_figure,
)

try:
    from adjustText import adjust_text

    _HAVE_ADJUST = True
except Exception:  # pragma: no cover
    _HAVE_ADJUST = False


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_table(paths, "s4", "significance_screen").copy()
    df = df[df["is_signed"].astype(bool)].copy()
    df["beta_signed"] = pd.to_numeric(df["beta_signed"], errors="coerce")
    df = df.dropna(subset=["beta_signed"]).reset_index(drop=True)
    df["chain"] = df["chain"].astype(str)
    df["_resid"] = df["canon_label"].map(lambda s: parse_canon_label(s)[1])

    fig, ax = plt.subplots(figsize=(DOUBLE_COL, 3.2))
    if df.empty:
        ax.text(0.5, 0.5, "no signed mechanisms", ha="center", va="center",
                transform=ax.transAxes, color=styles.NONSIG_COLOR)
        ax.set_axis_off()
        return save_figure(fig, paths, "Supplementary_Figure_S2")

    order = (df[["chain", "_resid", "canon_label"]].drop_duplicates()
             .sort_values(["chain", "_resid"]).reset_index(drop=True))
    xpos = {lab: i for i, lab in enumerate(order["canon_label"])}
    df["_x"] = df["canon_label"].map(xpos)

    for i, (chain, s, e) in enumerate(chain_bands(list(order["chain"]))):
        ax.axvspan(s - 0.5, e - 0.5, color=("#f4f4f4" if i % 2 == 0 else "#eaeef3"),
                   zorder=0, lw=0)
        ax.text((s + e - 1) / 2, 1.01, chain, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=styles.FS_LABEL, fontweight="bold")

    ax.axhline(0, color=styles.OKABE_ITO["black"], lw=0.7, zorder=1)
    sig = df["significant_fdr"].astype(bool).to_numpy()
    # colour by FDR significance (identical semantics to Figure 4); sign shown by y
    for label, mask, color in [("not significant", ~sig, styles.NONSIG_COLOR),
                               ("FDR-significant", sig, styles.SIG_COLOR)]:
        d = df[mask]
        ax.vlines(d["_x"], 0, d["beta_signed"], color=color, lw=0.6, alpha=0.6,
                  zorder=2)
        ax.scatter(d["_x"], d["beta_signed"], s=(28 if label.startswith("FDR") else 12),
                   c=color, alpha=0.9, edgecolors="white", linewidths=0.3, zorder=3,
                   label=label)
    cat = is_catalytic(df["domain"])
    cat_sig = cat & df["significant_fdr"].astype(bool)
    if cat_sig.any():
        ax.scatter(df["_x"][cat_sig], df["beta_signed"][cat_sig], s=44, facecolors="none",
                   edgecolors=styles.CATALYTIC_ACCENT, linewidths=0.75, alpha=0.6,
                   zorder=4,
                   label="catalytic + FDR")

    idx = objective_labels(df, score="beta_signed", catalytic=None,
                           significant=df["significant_fdr"].astype(bool))
    texts = [ax.text(df.loc[i, "_x"], df.loc[i, "beta_signed"],
                     f"{df.loc[i, 'serotype']} {df.loc[i, 'canon_label']}",
                     fontsize=styles.FS_ANNOT) for i in idx]
    if _HAVE_ADJUST and texts:
        adjust_text(texts, ax=ax, iter_lim=250,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#bdbdbd"))

    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_xlabel("canonical residue position (ordered by chain, number)")
    ax.set_ylabel("signed effect β (energy, a.u.)")
    if len(order) <= 30:
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order["canon_label"], rotation=90, fontsize=styles.FS_ANNOT)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.62), ncol=3,
              fontsize=styles.FS_ANNOT)
    ax.set_title("Signed dynamic effects are distributed along NS2B–NS3", pad=16, y=1.08)
    fig.subplots_adjust(top=0.82, bottom=0.34)
    return save_figure(fig, paths, "Supplementary_Figure_S2")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
