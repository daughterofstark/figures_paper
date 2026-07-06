"""Supplementary Figure S2 — Conservation versus dynamic importance.

Scatter of cross-serotype conservation (x) against signed effect β (y), one point
per signed (serotype, residue) mechanism. Point colour encodes domain, point size
encodes per-position reproducibility ρ, catalytic residues are ringed, and a
deterministic LOWESS trend is overlaid. Only the strongest outliers are labelled.
It answers: *are evolutionarily conserved residues dynamically important?*

Inputs (existing outputs, joined on ``canon_label``):
``outputs_s5/position_conservation.parquet`` (frac_reproducible, rho_residue_median,
is_catalytic_triad) and ``outputs_s4/significance_screen.parquet`` (beta_signed).
No quantities are recomputed; the LOWESS curve is a visual trend only.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import ONEHALF_COL, Paths
from utils import is_catalytic, load_table, lowess, save_figure

try:
    from adjustText import adjust_text

    _HAVE_ADJUST = True
except Exception:  # pragma: no cover
    _HAVE_ADJUST = False


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    cons = load_table(paths, "s5", "position_conservation")[
        ["canon_label", "frac_reproducible", "rho_residue_median", "is_catalytic_triad"]
    ].copy()
    sig = load_table(paths, "s4", "significance_screen").copy()
    sig = sig[sig["is_signed"].astype(bool)]
    df = sig.merge(cons, on="canon_label", how="inner")
    df["frac_reproducible"] = pd.to_numeric(df["frac_reproducible"], errors="coerce")
    df["beta_signed"] = pd.to_numeric(df["beta_signed"], errors="coerce")
    df["rho_residue_median"] = pd.to_numeric(df["rho_residue_median"], errors="coerce")
    df = df.dropna(subset=["frac_reproducible", "beta_signed"])

    fig, ax = plt.subplots(figsize=(ONEHALF_COL, ONEHALF_COL * 0.82))
    if df.empty:
        ax.text(0.5, 0.5, "no signed mechanisms", ha="center", va="center",
                transform=ax.transAxes, color=styles.OKABE_ITO["grey"])
        ax.set_axis_off()
        return save_figure(fig, paths, "Supplementary_Figure_S2")

    ax.axhline(0, color=styles.OKABE_ITO["grey"], lw=0.6, zorder=1)

    domains = list(dict.fromkeys(df["domain"].astype(str)))
    palette = list(styles.OKABE_ITO.values())
    dom_color = {d: palette[i % len(palette)] for i, d in enumerate(domains)}
    sizes = 12 + 60 * df["rho_residue_median"].fillna(0)

    # small deterministic x jitter so overlapping conservation values separate
    xj = df["frac_reproducible"] + (df.groupby("frac_reproducible").cumcount()
                                    - 0.0) * 0.002
    for d in domains:
        m = df["domain"].astype(str) == d
        ax.scatter(xj[m], df["beta_signed"][m], s=sizes[m],
                   c=[dom_color[d]], alpha=0.8, edgecolors="white",
                   linewidths=0.3, label=d, zorder=2)
    cat = is_catalytic(df["domain"]) | df["is_catalytic_triad"].astype(bool)
    if cat.any():
        ax.scatter(xj[cat], df["beta_signed"][cat], s=70, facecolors="none",
                   edgecolors=styles.OKABE_ITO["black"], linewidths=1.0, zorder=3)

    trend = lowess(df["frac_reproducible"].to_numpy(), df["beta_signed"].to_numpy())
    if trend is not None:
        ax.plot(trend[0], trend[1], color=styles.OKABE_ITO["black"], lw=1.5,
                ls=(0, (5, 2)), zorder=4, label="LOWESS trend")

    # annotate the strongest |β| outliers
    strongest = df.reindex(df["beta_signed"].abs().sort_values(ascending=False).index)
    texts = []
    for _, r in strongest.head(6).iterrows():
        texts.append(ax.text(r["frac_reproducible"], r["beta_signed"],
                             f"{r['serotype']} {r['canon_label']}",
                             fontsize=styles.FS_ANNOT))
    if _HAVE_ADJUST and texts:
        adjust_text(texts, ax=ax, iter_lim=200,
                    arrowprops=dict(arrowstyle="-", lw=0.4,
                                    color=styles.OKABE_ITO["grey"]))

    ax.set_xlabel("cross-serotype conservation (fraction reproducible)")
    ax.set_ylabel("signed effect β (energy, a.u.)")
    ax.set_title("Conservation vs dynamic importance")
    ax.margins(x=0.08, y=0.12)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=styles.FS_ANNOT, title="domain")
    return save_figure(fig, paths, "Supplementary_Figure_S2")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
