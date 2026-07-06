"""Figure 4 — Signed-effect volcano (redesign of S7 F4).

The original F4 forest plots hundreds of confidence intervals and is unreadable.
This replaces it with a **volcano-style summary** that scales to the full dataset:
signed effect size (β) on x, statistical strength (−log10 BH-adjusted p) on y,
FDR-significant points accented, catalytic residues ringed and labelled. No new
statistics are computed — β, p and the FDR flag are read verbatim from S4.

Consumes: ``outputs_s4/significance_screen.parquet`` (signed mechanisms + p-values).
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import ONEHALF_COL, Paths
from utils import is_catalytic, load_table, save_figure

try:
    from adjustText import adjust_text

    _HAVE_ADJUST = True
except Exception:  # pragma: no cover - optional dependency
    _HAVE_ADJUST = False


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_table(paths, "s4", "significance_screen").copy()
    df = df[df["is_signed"].astype(bool)].copy()
    df["beta_signed"] = pd.to_numeric(df["beta_signed"], errors="coerce")
    df["p_value_bh"] = pd.to_numeric(df["p_value_bh"], errors="coerce")
    df = df.dropna(subset=["beta_signed", "p_value_bh"])

    fig, ax = plt.subplots(figsize=(ONEHALF_COL, ONEHALF_COL * 0.72))

    if df.empty:
        ax.text(0.5, 0.5, "no signed mechanisms", ha="center", va="center",
                transform=ax.transAxes, color=styles.OKABE_ITO["grey"])
        ax.set_axis_off()
        return save_figure(fig, paths, "figure4_signed_effect_volcano")

    # clip p to avoid infinities; y = -log10(p_bh)
    pmin = float(df.loc[df["p_value_bh"] > 0, "p_value_bh"].min() or 1e-300)
    y = -np.log10(df["p_value_bh"].clip(lower=pmin / 10))
    x = df["beta_signed"].to_numpy()
    sig = df["significant_fdr"].astype(bool).to_numpy()
    cat = is_catalytic(df["domain"]).to_numpy()

    ax.axvline(0, color=styles.OKABE_ITO["grey"], lw=0.6, zorder=1)
    ax.axhline(-np.log10(0.05), color=styles.OKABE_ITO["grey"], lw=0.6,
               ls=(0, (4, 3)), zorder=1)

    ax.scatter(x[~sig], y[~sig], s=12, c=styles.NONSIG_COLOR, alpha=0.6,
               edgecolors="none", zorder=2, label="not FDR-significant")
    ax.scatter(x[sig], y[sig], s=16, c=styles.SIG_COLOR, alpha=0.85,
               edgecolors="none", zorder=3, label="FDR-significant")
    if cat.any():
        ax.scatter(x[cat], y[cat], s=40, facecolors="none",
                   edgecolors=styles.OKABE_ITO["black"], linewidths=1.0,
                   zorder=4, label="catalytic residue")

    # label catalytic + the strongest few effects
    df_lab = df.assign(_y=y.to_numpy())
    to_label = df_lab[cat]
    strongest = df_lab.reindex(df_lab["_y"].sort_values(ascending=False).index).head(6)
    label_rows = pd.concat([to_label, strongest]).drop_duplicates()
    texts = []
    for _, r in label_rows.iterrows():
        texts.append(ax.text(r["beta_signed"], r["_y"],
                             f"{r['serotype']} {r['canon_label']}",
                             fontsize=styles.FS_ANNOT))
    if _HAVE_ADJUST and texts:
        # adjustText defaults to a wall-clock stop (time_lim=1s), which makes the
        # result machine-speed-dependent. A fixed iter_lim makes it deterministic
        # (its internal RNG is already fixed).
        adjust_text(texts, ax=ax, iter_lim=200,
                    arrowprops=dict(arrowstyle="-", lw=0.4,
                                    color=styles.OKABE_ITO["grey"]))

    ax.set_xlabel("signed effect β (energy, a.u.)")
    ax.set_ylabel(r"$-\log_{10}$ (BH-adjusted $p$)")
    ax.set_title("Signed-effect volcano (coherent mechanisms)")
    ax.legend(loc="upper left", bbox_to_anchor=(0.0, 1.0))
    ax.margins(x=0.12, y=0.12)
    fig.tight_layout()
    return save_figure(fig, paths, "figure4_signed_effect_volcano")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
