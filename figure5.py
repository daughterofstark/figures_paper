"""Figure 5 — A conserved reproducible core is dynamically important.

Biological message : a core set of positions reproduces across all four serotypes,
    the serotypes cluster into a clear similarity structure at the dynamic level, and
    the most conserved positions carry the largest signed effects — conservation and
    dynamic importance coincide. This is the manuscript's central cross-serotype claim.
Why this figure exists : it is the paper's headline — it is the only figure that links
    reproducibility to evolutionary/functional importance across serotypes.
Generated from : outputs_s5/position_conservation.parquet (conservation per position),
    outputs_s7/F1_reproducibility_landscape.parquet (per-serotype ρ profiles → 4×4
    similarity), and outputs_s4/significance_screen.parquet (signed effects).
Placement : MAIN (headline). Merges the former conservation, serotype-matrix and
    conservation-vs-effect figures.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import CONSERVATION_ORDER, DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import (
    is_catalytic,
    load_s7_figure,
    load_table,
    lowess,
    objective_labels,
    panel_letter,
    parse_canon_label,
    save_figure,
)

try:
    from adjustText import adjust_text

    _HAVE_ADJUST = True
except Exception:  # pragma: no cover
    _HAVE_ADJUST = False


def _panel_conservation(ax, cons: pd.DataFrame) -> None:
    d = cons.copy()
    d["n_serotypes_reproducible"] = pd.to_numeric(d["n_serotypes_reproducible"],
                                                  errors="coerce")
    d["_key"] = d["canon_label"].map(parse_canon_label)
    d = d.sort_values(["n_serotypes_reproducible", "_key"],
                      ascending=[True, True]).reset_index(drop=True)
    y = np.arange(len(d))
    colors = [styles.CONSERVATION_COLORS.get(c, styles.NONSIG_COLOR)
              for c in d["conservation_class"]]
    ax.hlines(y, 0, d["n_serotypes_reproducible"], color="#dcdcdc", lw=0.5, zorder=1)
    ax.scatter(d["n_serotypes_reproducible"], y, s=22, c=colors, edgecolors="white",
               linewidths=0.4, zorder=2)
    cat = d[d.get("is_catalytic_triad", pd.Series(False, index=d.index)).astype(bool)]
    if not cat.empty:
        ax.scatter(cat["n_serotypes_reproducible"], cat.index, s=64,
                   facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                   linewidths=1.1, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(d["canon_label"], fontsize=styles.FS_ANNOT)
    ax.set_xlabel("serotypes reproducing (n)")
    ax.set_xlim(-0.2, max(4, float(d["n_serotypes_reproducible"].max() or 4)) + 0.3)
    ax.set_title("a conserved reproducible core", fontsize=styles.FS_LABEL)
    handles = [plt.Line2D([], [], marker="o", ls="none", mfc=styles.CONSERVATION_COLORS[c],
                          mec="white", ms=5, label=c.replace("reproducible_", ""))
               for c in CONSERVATION_ORDER if c in set(d["conservation_class"])]
    handles.append(plt.Line2D([], [], marker="o", ls="none", mfc="none",
                              mec=styles.CATALYTIC_ACCENT, ms=7, mew=1.1,
                              label="catalytic triad"))
    ax.legend(handles=handles, loc="lower right", fontsize=styles.FS_ANNOT)


def _panel_matrix(fig, cell, land: pd.DataFrame) -> None:
    land = land.copy()
    land["rho_residue"] = pd.to_numeric(land["rho_residue"], errors="coerce")
    mat = land.pivot_table(index="canon_label", columns="serotype",
                           values="rho_residue", aggfunc="mean")
    serotypes = [s for s in SEROTYPE_ORDER if s in mat.columns]
    serotypes += [s for s in mat.columns if s not in serotypes]
    mat = mat[serotypes]
    corr = mat.corr(method="pearson", min_periods=2)
    order = list(corr.columns)
    corr_vals = corr.to_numpy()

    axd = None
    if len(order) >= 3 and np.isfinite(corr_vals).all():
        dist = 1.0 - corr_vals
        np.fill_diagonal(dist, 0.0)
        dist = (dist + dist.T) / 2.0
        if float(np.nanmax(dist)) >= 1e-9:
            try:
                from scipy.cluster.hierarchy import dendrogram, linkage
                from scipy.spatial.distance import squareform

                z = linkage(squareform(dist, checks=False), method="average")
                sub = cell.subgridspec(2, 1, height_ratios=[1, 4], hspace=0.05)
                axd = fig.add_subplot(sub[0])
                dd = dendrogram(z, no_labels=True, color_threshold=0,
                                above_threshold_color="#9e9e9e", ax=axd)
                axd.set_axis_off()
                order = [corr.columns[i] for i in dd["leaves"]]
                ax = fig.add_subplot(sub[1])
            except Exception:
                ax = fig.add_subplot(cell)
        else:
            ax = fig.add_subplot(cell)
    else:
        ax = fig.add_subplot(cell)

    cm = corr.loc[order, order]
    im = ax.imshow(cm.to_numpy(float), cmap=styles.CORR_CMAP, vmin=-1, vmax=1,
                   aspect="equal")
    ax.set_xticks(range(len(order)))
    ax.set_yticks(range(len(order)))
    ax.set_xticklabels(order, rotation=45, ha="right", fontsize=styles.FS_ANNOT)
    ax.set_yticklabels(order, fontsize=styles.FS_ANNOT)
    for i in range(len(order)):
        for j in range(len(order)):
            v = cm.to_numpy()[i, j]
            if np.isfinite(v):
                ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                        fontsize=styles.FS_ANNOT,
                        color="white" if v < 0.35 else "black")
    # title on the top-most axis (dendrogram if present) so nothing overlaps
    top_ax = axd if axd is not None else ax
    top_ax.set_title("serotype similarity (ρ profiles)", fontsize=styles.FS_LABEL,
                     pad=4)
    cbar = fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    cbar.set_label("Pearson r", fontsize=styles.FS_TICK)
    cbar.set_ticks([-1, 0, 1])
    cbar.outline.set_linewidth(0.4)
    return axd if axd is not None else ax


def _panel_cons_vs_effect(ax, cons: pd.DataFrame, sig: pd.DataFrame) -> None:
    c = cons[["canon_label", "frac_reproducible", "rho_residue_median",
              "is_catalytic_triad"]].copy()
    s = sig[sig["is_signed"].astype(bool)].copy()
    df = s.merge(c, on="canon_label", how="inner")
    df["frac_reproducible"] = pd.to_numeric(df["frac_reproducible"], errors="coerce")
    df["beta_signed"] = pd.to_numeric(df["beta_signed"], errors="coerce")
    df = df.dropna(subset=["frac_reproducible", "beta_signed"]).reset_index(drop=True)
    if df.empty:
        ax.text(0.5, 0.5, "no signed mechanisms", ha="center", va="center",
                transform=ax.transAxes, color=styles.NONSIG_COLOR)
        return
    ax.axhline(0, color="#c8c8c8", lw=0.6, zorder=1)
    cat = is_catalytic(df["domain"]) | df["is_catalytic_triad"].astype(bool)
    xj = df["frac_reproducible"] + df.groupby("frac_reproducible").cumcount() * 0.003
    ax.scatter(xj[~cat], df["beta_signed"][~cat], s=16, c=styles.NONSIG_COLOR,
               alpha=0.8, edgecolors="white", linewidths=0.3, zorder=2,
               label="non-catalytic")
    ax.scatter(xj[cat], df["beta_signed"][cat], s=26, c=styles.CATALYTIC_ACCENT,
               alpha=0.9, edgecolors="white", linewidths=0.3, zorder=3,
               label="catalytic")
    tr = lowess(df["frac_reproducible"].to_numpy(), df["beta_signed"].to_numpy())
    if tr is not None:
        ax.plot(tr[0], tr[1], color=styles.OKABE_ITO["black"], lw=1.4,
                ls=(0, (5, 2)), zorder=4, label="trend (visual)")
    idx = objective_labels(df, score="beta_signed", catalytic=cat,
                           significant=df["significant_fdr"].astype(bool)
                           if "significant_fdr" in df else None)
    texts = [ax.text(df.loc[i, "frac_reproducible"], df.loc[i, "beta_signed"],
                     f"{df.loc[i, 'serotype']} {df.loc[i, 'canon_label']}",
                     fontsize=styles.FS_ANNOT) for i in idx]
    if _HAVE_ADJUST and texts:
        adjust_text(texts, ax=ax, iter_lim=200,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color="#bdbdbd"))
    ax.set_xlabel("cross-serotype conservation (fraction reproducible)")
    ax.set_ylabel("signed effect β")
    ax.margins(x=0.06, y=0.12)
    ax.set_title("conserved positions carry the larger effects", fontsize=styles.FS_LABEL)
    ax.legend(loc="upper left", fontsize=styles.FS_ANNOT, ncol=1)


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    cons = load_table(paths, "s5", "position_conservation")
    land = load_s7_figure(paths, "F1_reproducibility_landscape")
    sig = load_table(paths, "s4", "significance_screen")

    fig = plt.figure(figsize=(DOUBLE_COL, 5.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.05, 1.0], width_ratios=[1.35, 1.0],
                          hspace=0.42, wspace=0.34)
    axa = fig.add_subplot(gs[0, 0])
    _panel_conservation(axa, cons)
    panel_letter(axa, "A")

    top_b = _panel_matrix(fig, gs[0, 1], land)
    panel_letter(top_b, "B")

    axc = fig.add_subplot(gs[1, :])
    _panel_cons_vs_effect(axc, cons, sig)
    panel_letter(axc, "C")

    fig.suptitle(
        "A conserved reproducible core, shared across DENV1–4, is dynamically important",
        x=0.5, y=1.0)
    return save_figure(fig, paths, "figure5_cross_serotype")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
