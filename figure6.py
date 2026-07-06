"""Figure 6 — Variance composition (redesign of S7 F6).

Separates serotypes into small-multiple panels and groups domains biologically
(catalytic machinery first, then the rest), showing the τ² (between-replicate) vs
σ̄² (within-replicate sampling) split as 100%-stacked horizontal bars. The τ²
fraction is the reproducibility-limiting component, so a large orange bar flags a
replicate-dominated domain at a glance.

Consumes: ``outputs_s7/F6_variance_composition.parquet``.
"""
from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import CATALYTIC_DOMAINS, DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import load_s7_figure, panel_letter, save_figure

_LETTERS = "ABCDEFGH"


def _domain_order(domains: list[str]) -> list[str]:
    catalytic = [d for d in CATALYTIC_DOMAINS if d in domains]
    others = sorted(d for d in domains if d not in CATALYTIC_DOMAINS)
    return catalytic + others


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df = load_s7_figure(paths, "F6_variance_composition").copy()
    df["frac_tau2"] = pd.to_numeric(df["frac_tau2"], errors="coerce").fillna(0)
    df["frac_sigma2"] = pd.to_numeric(df["frac_sigma2"], errors="coerce").fillna(0)

    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))
    domains = _domain_order(sorted(set(df["domain"].astype(str))))
    n = max(1, len(serotypes))
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))

    fig, axes = plt.subplots(
        nrow, ncol, figsize=(DOUBLE_COL, 0.95 * nrow + 0.8),
        squeeze=False, sharex=True,
    )
    flat = axes.ravel()

    for k, sero in enumerate(serotypes):
        ax = flat[k]
        sub = df[df["serotype"] == sero].set_index("domain").reindex(domains)
        y = np.arange(len(domains))
        tau = sub["frac_tau2"].to_numpy()
        sig = sub["frac_sigma2"].to_numpy()
        ax.barh(y, tau, color=styles.VARIANCE_COLORS["tau2"], edgecolor="white",
                linewidth=0.4, label=r"$\tau^2$ (between-replicate)")
        ax.barh(y, sig, left=tau, color=styles.VARIANCE_COLORS["sigma2"],
                edgecolor="white", linewidth=0.4,
                label=r"$\bar{\sigma}^2$ (within-replicate)")
        ax.set_yticks(y)
        ylabels = [
            (f"★ {d}" if d in CATALYTIC_DOMAINS else d) for d in domains
        ]
        ax.set_yticklabels(ylabels, fontsize=styles.FS_ANNOT)
        for tick, d in zip(ax.get_yticklabels(), domains, strict=True):
            if d in CATALYTIC_DOMAINS:
                tick.set_color(styles.CATALYTIC_ACCENT)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_title(sero, fontsize=styles.FS_LABEL)
        panel_letter(ax, _LETTERS[k])

    for k in range(len(serotypes), len(flat)):
        flat[k].set_visible(False)
    for ax in flat[len(serotypes) - ncol: len(serotypes)]:
        ax.set_xlabel("variance fraction")

    handles, labels = flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", ncol=2,
               bbox_to_anchor=(0.5, 1.005))
    fig.suptitle("Domain variance composition (τ² vs σ̄²) by serotype",
                 y=1.05)
    fig.tight_layout()
    return save_figure(fig, paths, "figure6_variance_composition")


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
