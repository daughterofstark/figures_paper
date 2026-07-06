"""Figure 8 — Coherence vs reproducibility (redesign of S7 F8).

Improves the scatter's readability and labels the biologically important domains.
Points are domain × serotype summaries; coherent domains (directionally clean) are
accented, catalytic domains are ringed and labelled. Provisional ρ* and coherence
guides split the plane into the "reproducible & directionally clean" quadrant
(top-right) versus the rest — annotated as provisional (uncalibrated).

Consumes: ``outputs_s7/F8_coherence_vs_rho.parquet``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import pandas as pd

import styles
from figure_config import (
    CATALYTIC_DOMAINS,
    ONEHALF_COL,
    PROVISIONAL_COHERENCE_THRESHOLD,
    PROVISIONAL_RHO_STAR,
    Paths,
)
from utils import is_catalytic, load_s7_figure, save_figure

try:
    from adjustText import adjust_text

    _HAVE_ADJUST = True
except Exception:  # pragma: no cover
    _HAVE_ADJUST = False


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F8_coherence_vs_rho").copy()
    df["rho_domain"] = pd.to_numeric(df["rho_domain"], errors="coerce")
    df["coherence_domain"] = pd.to_numeric(df["coherence_domain"], errors="coerce")
    df = df.dropna(subset=["rho_domain", "coherence_domain"])

    fig, ax = plt.subplots(figsize=(ONEHALF_COL, ONEHALF_COL * 0.82))

    # provisional quadrant guides
    ax.axvline(PROVISIONAL_RHO_STAR, color=styles.OKABE_ITO["grey"], lw=0.6,
               ls=(0, (4, 3)), zorder=1)
    ax.axhline(PROVISIONAL_COHERENCE_THRESHOLD, color=styles.OKABE_ITO["grey"],
               lw=0.6, ls=(0, (4, 3)), zorder=1)
    ax.add_patch(plt.Rectangle(
        (PROVISIONAL_RHO_STAR, PROVISIONAL_COHERENCE_THRESHOLD),
        1 - PROVISIONAL_RHO_STAR, 1 - PROVISIONAL_COHERENCE_THRESHOLD,
        facecolor=styles.OKABE_ITO["green"], alpha=0.06, zorder=0))

    coherent = df["is_coherent"].astype(bool).to_numpy()
    for label, mask, color in [
        ("coherent", coherent, styles.COHERENT_COLOR),
        ("mixed", ~coherent, styles.MIXED_COLOR),
    ]:
        d = df[mask]
        ax.scatter(d["rho_domain"], d["coherence_domain"], s=22, c=color,
                   alpha=0.8, edgecolors="white", linewidths=0.4, label=label,
                   zorder=2)

    cat = df[is_catalytic(df["domain"])]
    if not cat.empty:
        ax.scatter(cat["rho_domain"], cat["coherence_domain"], s=52,
                   facecolors="none", edgecolors=styles.CATALYTIC_ACCENT,
                   linewidths=1.2, zorder=3, label="catalytic domain")

    # label catalytic domains (one label per unique domain, at its centroid)
    texts = []
    for dom, g in cat.groupby("domain"):
        texts.append(ax.text(g["rho_domain"].mean(), g["coherence_domain"].mean(),
                             str(dom), fontsize=styles.FS_ANNOT,
                             color=styles.CATALYTIC_ACCENT))
    if _HAVE_ADJUST and texts:
        adjust_text(texts, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.4,
                    color=styles.OKABE_ITO["grey"]))

    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xlabel("ρ (domain reproducibility)")
    ax.set_ylabel("coherence (directional cleanliness)")
    ax.set_title("Reproducibility vs coherence at the domain scale")
    ax.text(0.99, 0.02, "guides are provisional (uncalibrated)", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=styles.FS_ANNOT,
            color=styles.OKABE_ITO["grey"])
    ax.legend(loc="lower left", bbox_to_anchor=(0.0, 0.0))
    fig.tight_layout()
    return save_figure(fig, paths, "figure8_coherence_vs_rho")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
