"""Figure 6 — Reproducibility is limited by between-replicate variance.

Biological message : where reproducibility is imperfect, the limiting term is τ²
    (disagreement between replicates), not σ̄² (within-replicate sampling noise) —
    so the ceiling on reproducibility is set by replicate-to-replicate variability.
Why this figure exists : it explains *what limits* reproducibility, the mechanistic
    counterpart to the earlier "where/which/how much" figures.
Generated from : outputs_s7/F6_variance_composition.parquet (frac_tau2, frac_sigma2
    per serotype × (chain, domain)).
Placement : MAIN. Kept as the full τ²/σ̄² decomposition (both components shown) so the
    balance is directly interpretable per region and serotype.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import CATALYTIC_DOMAINS, DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import panel_letter, save_figure

_LETTERS = "ABCDEFGH"


def _region_order(df: pd.DataFrame) -> list[tuple[str, str]]:
    """Unique (chain, domain) regions, catalytic first. (chain, domain) is the
    S4 variance-budget key and is unique within a serotype, so a domain label that
    recurs across chains is kept distinct rather than collapsed."""
    regions = list(dict.fromkeys(
        zip(df["chain"].astype(str), df["domain"].astype(str), strict=True)))

    def key(cd: tuple[str, str]) -> tuple[int, str, str]:
        chain, dom = cd
        rank = (CATALYTIC_DOMAINS.index(dom) if dom in CATALYTIC_DOMAINS
                else len(CATALYTIC_DOMAINS))
        return (rank, chain, dom)

    return sorted(regions, key=key)


def _region_label(chain: str, domain: str) -> str:
    star = "* " if domain in CATALYTIC_DOMAINS else ""
    return f"{star}{domain} ({chain})"


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = pd.read_parquet(paths.s7("F6_variance_composition.parquet")).copy()
    df["chain"] = df["chain"].astype(str)
    df["domain"] = df["domain"].astype(str)
    df["frac_tau2"] = pd.to_numeric(df["frac_tau2"], errors="coerce").fillna(0)
    df["frac_sigma2"] = pd.to_numeric(df["frac_sigma2"], errors="coerce").fillna(0)

    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))
    regions = _region_order(df)
    n = max(1, len(serotypes))
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))

    fig, axes = plt.subplots(nrow, ncol, figsize=(DOUBLE_COL, 0.95 * nrow + 0.8),
                             squeeze=False, sharex=True)
    flat = axes.ravel()
    for k, sero in enumerate(serotypes):
        ax = flat[k]
        sub = (df[df["serotype"] == sero].set_index(["chain", "domain"])
               .reindex(regions))
        y = np.arange(len(regions))
        tau = sub["frac_tau2"].to_numpy()
        sig = sub["frac_sigma2"].to_numpy()
        ax.barh(y, tau, color=styles.VARIANCE_COLORS["tau2"], edgecolor="white",
                linewidth=0.4, label=r"$\tau^2$ (between-replicate)")
        ax.barh(y, sig, left=tau, color=styles.VARIANCE_COLORS["sigma2"],
                edgecolor="white", linewidth=0.4,
                label=r"$\bar{\sigma}^2$ (within-replicate)")
        ax.set_yticks(y)
        ax.set_yticklabels([_region_label(c, d) for c, d in regions],
                           fontsize=styles.FS_ANNOT)
        for tick, (_c, d) in zip(ax.get_yticklabels(), regions, strict=True):
            if d in CATALYTIC_DOMAINS:
                tick.set_color(styles.CATALYTIC_ACCENT)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_title(sero, fontsize=styles.FS_LABEL)
        panel_letter(ax, _LETTERS[k])

    for k in range(len(serotypes), len(flat)):
        flat[k].set_visible(False)
    for ax in flat[max(0, len(serotypes) - ncol):len(serotypes)]:
        ax.set_xlabel("variance fraction")

    handles, labels = flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 1.0),
               fontsize=styles.FS_ANNOT)
    fig.suptitle("Reproducibility is limited by between-replicate variance (τ²)",
                 y=1.05)
    fig.tight_layout()
    return save_figure(fig, paths, "figure6_variance_composition")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
