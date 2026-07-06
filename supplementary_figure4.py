"""Supplementary Figure S4 — Dynamic importance along the sequence.

A Manhattan / lollipop plot of the signed effect β at each canonical residue: x is
the canonical residue position, y is β, colour encodes sign (increase vs decrease),
and point size encodes FDR significance. A zero line is drawn, catalytic residues
are ringed, and only the top-10 positive and top-10 negative effects are labelled
(overlap-avoided). This replaces the unreadable forest-plot concept.

Input: ``outputs_s4/significance_screen.parquet`` (beta_signed, significant_fdr,
is_signed, chain, domain). Catalytic residues are flagged from the existing domain
annotation. Nothing is recomputed.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import DOUBLE_COL, Paths
from utils import chain_bands, is_catalytic, load_table, panel_letter, parse_canon_label, save_figure

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
    df = df.dropna(subset=["beta_signed"])
    df["chain"] = df["chain"].astype(str)
    df["_resid"] = df["canon_label"].map(lambda s: parse_canon_label(s)[1])

    fig, ax = plt.subplots(figsize=(DOUBLE_COL, 3.2))
    if df.empty:
        ax.text(0.5, 0.5, "no signed mechanisms", ha="center", va="center",
                transform=ax.transAxes, color=styles.OKABE_ITO["grey"])
        ax.set_axis_off()
        return save_figure(fig, paths, "Supplementary_Figure_S4")

    # one column per canonical residue, ordered by (chain, resid); serotypes stack
    order = (
        df[["chain", "_resid", "canon_label"]].drop_duplicates()
        .sort_values(["chain", "_resid"]).reset_index(drop=True)
    )
    xpos = {lab: i for i, lab in enumerate(order["canon_label"])}
    df["_x"] = df["canon_label"].map(xpos)

    # chain background bands
    bands = chain_bands(list(order["chain"]))
    for i, (chain, s, e) in enumerate(bands):
        ax.axvspan(s - 0.5, e - 0.5, color=("#f2f2f2" if i % 2 == 0 else "#e9eef4"),
                   zorder=0, lw=0)
        ax.text((s + e - 1) / 2, 1.01, chain, transform=ax.get_xaxis_transform(),
                ha="center", va="bottom", fontsize=styles.FS_LABEL,
                fontweight="bold")

    ax.axhline(0, color=styles.OKABE_ITO["black"], lw=0.7, zorder=1)
    pos = df["beta_signed"] >= 0
    sig = df["significant_fdr"].astype(bool)
    sizes = np.where(sig, 34, 12)
    for label, mask, color in [
        ("increase (β > 0)", pos, styles.SIG_COLOR),
        ("decrease (β < 0)", ~pos, styles.COHERENT_COLOR),
    ]:
        d = df[mask]
        ax.vlines(d["_x"], 0, d["beta_signed"], color=color, lw=0.6, alpha=0.6,
                  zorder=2)
        ax.scatter(d["_x"], d["beta_signed"], s=sizes[mask.to_numpy()], c=color,
                   alpha=0.9, edgecolors="white", linewidths=0.3, zorder=3,
                   label=label)

    cat = is_catalytic(df["domain"])
    if cat.any():
        ax.scatter(df["_x"][cat], df["beta_signed"][cat], s=70, facecolors="none",
                   edgecolors=styles.OKABE_ITO["black"], linewidths=1.0, zorder=4,
                   label="catalytic residue")

    # annotate top-10 positive and top-10 negative
    top_pos = df[pos].sort_values("beta_signed", ascending=False).head(10)
    top_neg = df[~pos].sort_values("beta_signed", ascending=True).head(10)
    texts = []
    for _, r in pd.concat([top_pos, top_neg]).iterrows():
        texts.append(ax.text(r["_x"], r["beta_signed"],
                             f"{r['serotype']} {r['canon_label']}",
                             fontsize=styles.FS_ANNOT))
    if _HAVE_ADJUST and texts:
        adjust_text(texts, ax=ax, iter_lim=250,
                    arrowprops=dict(arrowstyle="-", lw=0.4,
                                    color=styles.OKABE_ITO["grey"]))

    ax.set_xlim(-0.6, len(order) - 0.4)
    ax.set_xlabel("canonical residue position (ordered by chain, number)")
    ax.set_ylabel("signed effect β (energy, a.u.)")
    ax.set_title("Dynamic importance along the sequence", pad=16, y=1.08)
    if len(order) <= 30:
        ax.set_xticks(range(len(order)))
        ax.set_xticklabels(order["canon_label"], rotation=90,
                           fontsize=styles.FS_ANNOT)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.62), ncol=3,
              fontsize=styles.FS_ANNOT)
    fig.subplots_adjust(top=0.82, bottom=0.34)
    return save_figure(fig, paths, "Supplementary_Figure_S4")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
