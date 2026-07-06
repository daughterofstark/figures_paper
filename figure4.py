"""Figure 4 — Significant signed effects concentrate in catalytic regions.

Biological message : only a minority of signed dynamic effects survive FDR control,
    and those that do are enriched in the catalytic machinery.
Why this figure exists : it summarises effect magnitude and significance together
    (replacing an unreadable forest of hundreds of CIs) and shows the functional
    localisation of the significant effects.
Generated from : outputs_s4/significance_screen.parquet (beta_signed, p_value_bh,
    significant_fdr; catalytic flag from the domain annotation).
Placement : MAIN.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import ONEHALF_COL, Paths
from utils import is_catalytic, load_table, objective_labels, save_figure

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
    df["p_value_bh"] = pd.to_numeric(df["p_value_bh"], errors="coerce")
    df = df.dropna(subset=["beta_signed", "p_value_bh"]).reset_index(drop=True)

    fig, ax = plt.subplots(figsize=(ONEHALF_COL, ONEHALF_COL * 0.74))
    if df.empty:
        ax.text(0.5, 0.5, "no signed mechanisms", ha="center", va="center",
                transform=ax.transAxes, color=styles.NONSIG_COLOR)
        ax.set_axis_off()
        return save_figure(fig, paths, "figure4_signed_effect_volcano")

    pmin = float(df.loc[df["p_value_bh"] > 0, "p_value_bh"].min() or 1e-300)
    y = -np.log10(df["p_value_bh"].clip(lower=pmin / 10))
    x = df["beta_signed"].to_numpy()
    sig = df["significant_fdr"].astype(bool).to_numpy()
    cat = is_catalytic(df["domain"])

    ax.axvline(0, color="#c8c8c8", lw=0.6, zorder=1)
    ax.axhline(-np.log10(0.05), color="#c8c8c8", lw=0.6, ls=(0, (4, 3)), zorder=1)
    ax.text(ax.get_xlim()[0], -np.log10(0.05), " FDR 0.05", va="bottom", ha="left",
            fontsize=styles.FS_ANNOT, color="#9e9e9e")

    ax.scatter(x[~sig], y[~sig], s=12, c=styles.NONSIG_COLOR, alpha=0.7,
               edgecolors="none", zorder=2, label="not significant")
    ax.scatter(x[sig], y[sig], s=17, c=styles.SIG_COLOR, alpha=0.9,
               edgecolors="none", zorder=3, label="FDR-significant")
    if cat.any():
        ax.scatter(x[cat.to_numpy()], y[cat.to_numpy()], s=42, facecolors="none",
                   edgecolors=styles.CATALYTIC_ACCENT, linewidths=1.1, zorder=4,
                   label="catalytic")

    # single global label rule: catalytic OR (FDR-significant & top-decile |β|)
    df_lab = df.assign(_y=y.to_numpy())
    idx = objective_labels(df_lab, score="beta_signed", catalytic=cat,
                           significant=df_lab["significant_fdr"].astype(bool))
    texts = [ax.text(r["beta_signed"], r["_y"], f"{r['serotype']} {r['canon_label']}",
                     fontsize=styles.FS_ANNOT) for _, r in df_lab.loc[idx].iterrows()]
    if _HAVE_ADJUST and texts:
        adjust_text(texts, ax=ax, iter_lim=200,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#bdbdbd"))

    ax.set_xlabel("signed effect β (energy, a.u.)")
    ax.set_ylabel(r"$-\log_{10}$ (BH-adjusted $p$)")
    ax.margins(x=0.12, y=0.12)
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=styles.FS_ANNOT)
    fig.suptitle("Significant signed effects concentrate in catalytic regions",
                 x=0.5, y=1.0)
    return save_figure(fig, paths, "figure4_signed_effect_volcano")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
