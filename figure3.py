"""Figure 3 — Catalytic domains are coherent but not the highest-ρ regions.

Biological message : domain-scale ρ is highest in broad non-catalytic or unassigned
    regions, whereas catalytic domains are moderate in ρ but comparatively coherent
    in signed direction.
Why this figure exists : it identifies *which structural module* carries the signal,
    combining the cross-serotype ρ landscape with the ρ–coherence relationship so the
    two are read as one result rather than two.
Generated from : outputs_s7/F3_domain_serotype_rho_heatmap.parquet (ρ per
    serotype × domain) and outputs_s7/F8_coherence_vs_rho.parquet (coherence vs ρ).
Placement : MAIN. Merges the former domain heatmap and coherence-vs-ρ figures.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

import styles
from figure_config import (
    CATALYTIC_DOMAINS,
    DOUBLE_COL,
    PROVISIONAL_COHERENCE_THRESHOLD,
    PROVISIONAL_RHO_STAR,
    SEROTYPE_ORDER,
    Paths,
)
from utils import is_catalytic, load_s7_figure, panel_letter, save_figure


def _text_color(cmap, norm, v: float) -> str:
    r, g, b, _ = cmap(norm(v))
    return "white" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.55 else "black"


def _panel_heatmap(ax, df: pd.DataFrame, serotypes: list[str]):
    df = df.copy()
    df["chain"] = df["chain"].astype(str)
    df["domain"] = df["domain"].astype(str)
    df["rho_domain"] = pd.to_numeric(df["rho_domain"], errors="coerce")
    mat = df.pivot_table(index=["chain", "domain"], columns="serotype",
                         values="rho_domain", aggfunc="mean").reindex(columns=serotypes)
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
    regions = list(mat.index)
    norm = Normalize(0, 1)
    cmap = plt.get_cmap(styles.SEQ_CMAP)
    im = ax.imshow(mat.to_numpy(float), aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(serotypes)))
    ax.set_xticklabels(serotypes, rotation=0)
    ax.set_yticks(range(len(regions)))
    ax.set_yticklabels(
        [(f"{d} ({c})") for c, d in regions], fontsize=styles.FS_ANNOT)
    for tick, (_c, d) in zip(ax.get_yticklabels(), regions, strict=True):
        if d in CATALYTIC_DOMAINS:
            tick.set_color(styles.CATALYTIC_ACCENT)
            tick.set_fontweight("bold")
    vals = mat.to_numpy(float)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if np.isfinite(vals[i, j]):
                ax.text(j, i, f"{vals[i, j]:.2f}", ha="center", va="center",
                        fontsize=styles.FS_ANNOT,
                        color=_text_color(cmap, norm, vals[i, j]))
    ax.set_title("domain-scale ρ is heterogeneous", fontsize=styles.FS_LABEL)
    return im, norm, cmap


def _panel_scatter(ax, df: pd.DataFrame, serotypes: list[str]):
    df = df.copy()
    df["rho_domain"] = pd.to_numeric(df["rho_domain"], errors="coerce")
    df["coherence_domain"] = pd.to_numeric(df["coherence_domain"], errors="coerce")
    df = df.dropna(subset=["rho_domain", "coherence_domain"])
    # provisional guides — faint, explicitly uncalibrated
    ax.axvline(PROVISIONAL_RHO_STAR, color="#d9d9d9", lw=0.6, ls=(0, (4, 3)), zorder=1)
    ax.axhline(PROVISIONAL_COHERENCE_THRESHOLD, color="#d9d9d9", lw=0.6,
               ls=(0, (4, 3)), zorder=1)
    cat_mask = is_catalytic(df["domain"])
    non = df[~cat_mask]
    ax.scatter(non["rho_domain"], non["coherence_domain"], s=22,
               c="#bdbdbd", alpha=0.62, edgecolors="white", linewidths=0.3,
               label="other domains", zorder=2)
    cat = df[cat_mask]
    if not cat.empty:
        ax.scatter(cat["rho_domain"], cat["coherence_domain"], s=50,
                   c=styles.CATALYTIC_ACCENT, alpha=0.88, edgecolors="white",
                   linewidths=0.5, zorder=3, label="catalytic domains")
        for dom, g in cat.groupby("domain", sort=False):
            x = pd.to_numeric(g["rho_domain"], errors="coerce").median()
            y = pd.to_numeric(g["coherence_domain"], errors="coerce").median()
            ax.text(x, y + 0.045, dom.replace("Catalytic ", "").replace(" Loop", ""),
                    ha="center", va="bottom", fontsize=styles.FS_ANNOT,
                    color=styles.CATALYTIC_ACCENT)
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("reproducibility ρ (domain)")
    ax.set_ylabel("coherence (directional cleanliness)")
    ax.set_title("coherence distinguishes catalytic domains", fontsize=styles.FS_LABEL)
    ax.text(0.99, 0.02, "guides provisional (uncalibrated)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=styles.FS_ANNOT, color="#9e9e9e")


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    heat = load_s7_figure(paths, "F3_domain_serotype_rho_heatmap")
    scat = load_s7_figure(paths, "F8_coherence_vs_rho")
    serotypes = [s for s in SEROTYPE_ORDER if s in set(heat["serotype"])]
    serotypes += sorted(set(heat["serotype"].astype(str)) - set(serotypes))

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 2.9),
        gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.32})
    im, norm, cmap = _panel_heatmap(axa, heat, serotypes)
    panel_letter(axa, "A")
    cbar = fig.colorbar(im, ax=axa, fraction=0.045, pad=0.03)
    cbar.set_label("ρ", fontsize=styles.FS_TICK)
    cbar.set_ticks([0, 0.5, 1.0])
    cbar.outline.set_linewidth(0.4)

    _panel_scatter(axb, scat, serotypes)
    axb.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0), fontsize=styles.FS_ANNOT,
               ncol=1)
    panel_letter(axb, "B", dx=-0.06)

    fig.suptitle("Catalytic domains are coherent but not the highest-ρ regions",
                 x=0.5, y=1.02)
    return save_figure(fig, paths, "figure3_catalytic_core")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
