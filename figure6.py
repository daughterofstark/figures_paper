"""Figure 6: fitted variance composition under the implemented K=3 model."""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import styles
from figure_config import DOUBLE_COL, SEROTYPE_ORDER, Paths
from phase3j_support import (
    group_support_counts,
    load_configured_group_metadata,
    sha256,
    write_json,
)
from utils import load_s7_figure, panel_letter, save_figure

_LETTERS = "ABCDEFGH"
TITLE = "Fitted variance composition under the K=3 diagonal-covariance model"
MODEL_NOTE = (
    "Propagated within-run uncertainty assumes zero item covariance; "
    "variance components are weakly identified at K=3."
)
FIGURE_STRINGS = (
    TITLE,
    MODEL_NOTE,
    "Fitted between-run component",
    "Propagated within-run component",
)


def _region_order(df: pd.DataFrame) -> list[tuple[str, str]]:
    """Deterministic configured-group ordering, flagged groups first."""
    regions = list(dict.fromkeys(
        zip(df["chain"].astype(str), df["domain"].astype(str), strict=True)))

    def key(cd: tuple[str, str]) -> tuple[int, str, str]:
        chain, dom = cd
        rows = df[(df["chain"] == chain) & (df["domain"] == dom)]
        rank = 0 if rows["is_configured_catalytic_group"].all() else 1
        return (rank, chain, dom)

    return sorted(regions, key=key)


def prepare_data(paths: Paths) -> tuple[pd.DataFrame, dict[str, str]]:
    """Join frozen variance rows to exact configured-group support metadata."""
    df = load_s7_figure(paths, "F6_variance_composition").copy()
    required = {
        "serotype", "chain", "domain", "frac_tau2", "frac_sigma2", "variance_regime"
    }
    if missing := required - set(df.columns):
        raise ValueError(f"Figure 6 input missing columns: {sorted(missing)}")
    keys = ["serotype", "chain", "domain"]
    if df.duplicated(keys).any():
        raise ValueError("Figure 6 input has duplicate configured-group rows")
    for col in ("frac_tau2", "frac_sigma2"):
        df[col] = pd.to_numeric(df[col], errors="raise")
    if df[["frac_tau2", "frac_sigma2"]].isna().any().any():
        raise ValueError("Figure 6 variance fractions contain non-finite values")
    if not np.allclose(df["frac_tau2"] + df["frac_sigma2"], 1.0, atol=1e-8):
        raise ValueError("Figure 6 variance fractions do not sum to one")
    item = load_s7_figure(paths, "F1_reproducibility_landscape")
    counts = group_support_counts(item)
    metadata, provenance = load_configured_group_metadata(counts)
    keep = keys + ["support_size"]
    if "is_configured_catalytic_group" in metadata:
        keep.append("is_configured_catalytic_group")
    df = df.merge(metadata[keep], on=keys, how="left", validate="one_to_one")
    if df["support_size"].isna().any():
        raise ValueError("Figure 6 has configured groups without support-size metadata")
    df["support_size"] = df["support_size"].astype(int)
    if "is_configured_catalytic_group" not in df:
        df["is_configured_catalytic_group"] = False
    df["is_configured_catalytic_group"] = (
        df["is_configured_catalytic_group"].fillna(False).astype(bool)
    )
    return df, provenance


def _region_label(chain: str, group: str, support_size: int) -> str:
    return f"{group} ({chain}; n={support_size})"


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    df, metadata_provenance = prepare_data(paths)
    df["chain"] = df["chain"].astype(str)
    df["domain"] = df["domain"].astype(str)

    serotypes = [s for s in SEROTYPE_ORDER if s in set(df["serotype"])]
    serotypes += sorted(set(df["serotype"].astype(str)) - set(serotypes))
    regions = _region_order(df)
    n = max(1, len(serotypes))
    ncol = 2 if n > 1 else 1
    nrow = int(np.ceil(n / ncol))

    fig, axes = plt.subplots(nrow, ncol, figsize=(DOUBLE_COL, 1.25 * nrow + 1.0),
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
        labels = [
            _region_label(c, d, int(sub.loc[(c, d), "support_size"]))
            for c, d in regions
        ]
        ax.set_yticklabels(labels,
                           fontsize=styles.FS_ANNOT)
        for tick, (_c, d) in zip(ax.get_yticklabels(), regions, strict=True):
            if bool(sub.loc[(_c, d), "is_configured_catalytic_group"]):
                tick.set_color(styles.CATALYTIC_ACCENT)
        ax.invert_yaxis()
        ax.set_xlim(0, 1)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_title(sero, loc="left", fontsize=styles.FS_LABEL)
        panel_letter(ax, _LETTERS[k])

    for k in range(len(serotypes), len(flat)):
        flat[k].set_visible(False)
    for ax in flat[max(0, len(serotypes) - ncol):len(serotypes)]:
        ax.set_xlabel("Fitted variance fraction")

    handles, labels = flat[0].get_legend_handles_labels()
    labels = ["Fitted between-run component", "Propagated within-run component"]
    fig.legend(handles, labels, loc="upper center", ncol=2, bbox_to_anchor=(0.5, 0.94),
               fontsize=styles.FS_ANNOT)
    fig.suptitle(TITLE, y=0.995)
    fig.text(0.5, 0.955, MODEL_NOTE, ha="center", va="top",
             fontsize=styles.FS_ANNOT, color="#5f5f5f")
    fig.tight_layout(rect=(0, 0, 1, 0.90))
    written = save_figure(fig, paths, "figure6_variance_composition")
    regime_counts = {
        str(k): int(v)
        for k, v in df["variance_regime"].value_counts().sort_index().items()
    }
    rows = []
    for row in df.sort_values(["serotype", "chain", "domain"]).itertuples(index=False):
        rows.append(
            {
                "dataset": row.serotype,
                "chain": row.chain,
                "configured_group": row.domain,
                "support_size": int(row.support_size),
                "frac_tau2": float(row.frac_tau2),
                "frac_sigma2": float(row.frac_sigma2),
                "variance_regime": str(row.variance_regime),
            }
        )
    write_json(
        paths.output / "figure6_numeric_audit.json",
        {
            "K": 3,
            "classification_rule": {
                "replicate_dominated": "frac_tau2 >= 0.6",
                "sampling_dominated": "frac_tau2 <= 0.4",
                "balanced": "0.4 < frac_tau2 < 0.6",
            },
            "regime_counts": regime_counts,
            "metadata": metadata_provenance,
            "rows": rows,
            "source": str(paths.s7("F6_variance_composition.parquet")),
            "source_sha256": sha256(paths.s7("F6_variance_composition.parquet")),
            "zero_covariance_approximation": True,
            "weak_identification_at_K3": True,
        },
    )
    (paths.output / "figure6_caption_draft.md").write_text(
        "**Figure 6. Fitted variance composition under the K=3 "
        "diagonal-covariance model.** Bars partition each configured-group fit "
        "into its fitted between-run component and propagated within-run "
        "component; labels give exact support size n. Regimes are bookkeeping "
        "classes defined by the between-run fraction: at least 0.6, "
        "replicate-dominated; at most 0.4, sampling-dominated; otherwise "
        "balanced. With K=3, heterogeneity and the fitted variance components "
        "are weakly identified. Propagated group uncertainty assumes zero "
        "covariance among item estimates, so the composition is model-based "
        "and does not prescribe a causal remedy.\n"
    )
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
