"""Figure 3: configured-group magnitude rho and profile coherence.

Panel A shows frozen configured-group magnitude rho with exact support sizes.
Panel B plots magnitude rho against profile coherence: the mean over runs of the
within-run sign-homogeneity ratio |sum z| / sum |z| (the production ``coherence``
field written by ``statistics/reproducibility.py`` and carried into the profile by
``replicate/aggregator.build_profiles``). Profile coherence is NOT the pooled
coherence used for signed-analysis eligibility (``mech_coherence``, computed in
``replicate/aggregator._signed_fit``). The 0.6 line reproduces the pooled
eligibility value as a visual guide only.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.colors import Normalize

import styles
from figure_config import (
    DOUBLE_COL,
    PROVISIONAL_COHERENCE_THRESHOLD,
    PROVISIONAL_RHO_STAR,
    SEROTYPE_ORDER,
    Paths,
)
from phase3j_support import (
    group_support_counts,
    load_configured_group_metadata,
    sha256,
    write_json,
)
from utils import load_s7_figure, panel_letter, save_figure

TITLE = "Configured-group magnitude stability and profile coherence differ"
PANEL_TITLES = (
    "Magnitude rho by configured group",
    "Profile coherence versus magnitude rho",
)
Y_LABEL = "Profile coherence (mean within-run sign homogeneity)"
REFERENCE_NOTE = (
    "Guides: rho 0.5 provisional reference; "
    "0.6 = pooled-coherence eligibility value (visual guide only)"
)
FIGURE_STRINGS = (
    TITLE,
    *PANEL_TITLES,
    "Configured-group magnitude rho",
    Y_LABEL,
    "provisional rho reference",
    "pooled-coherence eligibility value",
)
COHERENCE_STATISTIC = "profile coherence: mean over runs of |sum z| / sum |z|"


def _text_color(cmap, norm, v: float) -> str:
    r, g, b, _ = cmap(norm(v))
    return "white" if (0.299 * r + 0.587 * g + 0.114 * b) < 0.55 else "black"


def prepare_data(paths: Paths) -> tuple[pd.DataFrame, dict[str, str]]:
    """Join frozen Figure 3 rows to exact configured-group support sizes."""
    heat = load_s7_figure(paths, "F3_domain_serotype_rho_heatmap").copy()
    scat = load_s7_figure(paths, "F8_coherence_vs_rho").copy()
    item = load_s7_figure(paths, "F1_reproducibility_landscape").copy()
    keys = ["serotype", "chain", "domain"]
    for name, frame, value_cols in (
        ("heatmap", heat, {"rho_domain"}),
        ("coherence", scat, {"rho_domain", "coherence_domain"}),
    ):
        required = set(keys) | value_cols
        if missing := required - set(frame.columns):
            raise ValueError(f"Figure 3 {name} input missing columns: {sorted(missing)}")
        if frame.duplicated(keys).any():
            raise ValueError(f"Figure 3 {name} input has duplicate configured-group rows")
    merged = heat.merge(
        scat[keys + ["rho_domain", "coherence_domain", "is_coherent"]],
        on=keys,
        how="outer",
        validate="one_to_one",
        suffixes=("_heat", "_scatter"),
        indicator=True,
    )
    if not (merged["_merge"] == "both").all():
        raise ValueError("Figure 3 heatmap/coherence configured-group keys differ")
    if not np.allclose(
        pd.to_numeric(merged["rho_domain_heat"], errors="raise"),
        pd.to_numeric(merged["rho_domain_scatter"], errors="raise"),
        rtol=0,
        atol=1e-12,
    ):
        raise ValueError("Figure 3 rho values differ between frozen source tables")
    merged = merged.drop(columns=["_merge", "rho_domain_scatter"]).rename(
        columns={"rho_domain_heat": "rho_domain"}
    )
    counts = group_support_counts(item)
    metadata, provenance = load_configured_group_metadata(counts)
    keep = keys + ["support_size"]
    if "is_configured_catalytic_group" in metadata:
        keep.append("is_configured_catalytic_group")
    merged = merged.merge(metadata[keep], on=keys, how="left", validate="one_to_one")
    if merged["support_size"].isna().any():
        raise ValueError("Figure 3 has configured groups without support-size metadata")
    merged["support_size"] = merged["support_size"].astype(int)
    if "is_configured_catalytic_group" not in merged:
        merged["is_configured_catalytic_group"] = False
    merged["is_configured_catalytic_group"] = (
        merged["is_configured_catalytic_group"].fillna(False).astype(bool)
    )
    return merged, provenance


def _panel_heatmap(ax, df: pd.DataFrame, serotypes: list[str]):
    df = df.copy()
    df["chain"] = df["chain"].astype(str)
    df["domain"] = df["domain"].astype(str)
    df["rho_domain"] = pd.to_numeric(df["rho_domain"], errors="coerce")
    mat = df.pivot(index=["chain", "domain"], columns="serotype",
                   values="rho_domain").reindex(columns=serotypes)
    mat = mat.loc[mat.mean(axis=1).sort_values(ascending=False).index]
    regions = list(mat.index)
    norm = Normalize(0, 1)
    cmap = plt.get_cmap(styles.SEQ_CMAP)
    im = ax.imshow(mat.to_numpy(float), aspect="auto", cmap=cmap, norm=norm)
    ax.set_xticks(range(len(serotypes)))
    ax.set_xticklabels(serotypes, rotation=0)
    ax.set_yticks(range(len(regions)))
    support = df.set_index(["chain", "domain", "serotype"])["support_size"]
    labels = [f"{group} ({chain})" for chain, group in regions]
    ax.set_yticklabels(labels, fontsize=styles.FS_ANNOT)
    for tick, (chain, group) in zip(ax.get_yticklabels(), regions, strict=True):
        group_rows = df[(df["chain"] == chain) & (df["domain"] == group)]
        if group_rows["is_configured_catalytic_group"].all():
            tick.set_color(styles.CATALYTIC_ACCENT)
            tick.set_fontweight("bold")
    vals = mat.to_numpy(float)
    for i in range(vals.shape[0]):
        for j in range(vals.shape[1]):
            if np.isfinite(vals[i, j]):
                size = int(support.loc[(regions[i][0], regions[i][1], serotypes[j])])
                # three decimals so values either side of the 0.5 reference
                # (e.g. 0.495 vs 0.502) are not both displayed as "0.50"
                ax.text(j, i, f"{vals[i, j]:.3f}\nn={size}",
                        ha="center", va="center",
                        fontsize=styles.FS_ANNOT,
                        color=_text_color(cmap, norm, vals[i, j]))
    ax.set_title(PANEL_TITLES[0], loc="left", fontsize=styles.FS_LABEL)
    return im, norm, cmap


def _panel_scatter(ax, df: pd.DataFrame):
    df = df.copy()
    df["rho_domain"] = pd.to_numeric(df["rho_domain"], errors="coerce")
    df["coherence_domain"] = pd.to_numeric(df["coherence_domain"], errors="coerce")
    df = df.dropna(subset=["rho_domain", "coherence_domain"])
    # visual guides only: rho 0.5 is provisional; 0.6 reproduces the pooled-
    # coherence eligibility value, which is a different statistic from the
    # profile coherence plotted on this axis
    ax.axvline(PROVISIONAL_RHO_STAR, color="#d9d9d9", lw=0.6, ls=(0, (4, 3)), zorder=1)
    ax.axhline(PROVISIONAL_COHERENCE_THRESHOLD, color="#d9d9d9", lw=0.6,
               ls=(0, (4, 3)), zorder=1)
    cat_mask = df["is_configured_catalytic_group"].astype(bool)
    non = df[~cat_mask]
    sizes = 14 + 1.8 * np.sqrt(non["support_size"].to_numpy(float))
    ax.scatter(non["rho_domain"], non["coherence_domain"], s=sizes,
               c="#737373", alpha=0.58, edgecolors="white", linewidths=0.3,
               label="other configured groups", zorder=2)
    cat = df[cat_mask]
    if not cat.empty:
        cat_sizes = 17 + 1.8 * np.sqrt(cat["support_size"].to_numpy(float))
        ax.scatter(cat["rho_domain"], cat["coherence_domain"], s=cat_sizes,
                   facecolors="white", edgecolors=styles.CATALYTIC_ACCENT,
                   linewidths=0.9, zorder=3, label="configured catalytic groups")
    ax.set_xlim(0, 1.02)
    ax.set_ylim(0, 1.02)
    ax.set_xticks([0, 0.5, 1.0])
    ax.set_yticks([0, 0.5, 1.0])
    ax.set_xlabel("Configured-group magnitude rho")
    ax.set_ylabel(Y_LABEL)
    ax.set_title(PANEL_TITLES[1], loc="left", fontsize=styles.FS_LABEL)


def _caption_text() -> str:
    return (
        "**Figure 3. Configured-group magnitude stability and profile coherence "
        "differ.** (A) Frozen magnitude rho for each of the 32 dataset-by-"
        "configured-group rows; each cell gives rho and the exact support size n. "
        "(B) Magnitude rho is shown against the implemented profile-coherence "
        "statistic, defined as the mean of the within-run sign-homogeneity values. "
        "Point area varies modestly with support size, and outlined points identify "
        "configured catalytic groups. The rho = 0.5 line is the provisional "
        "reference. The 0.6 line reproduces the value of the pooled-coherence "
        "eligibility reference and is drawn only as a visual guide. Signed-analysis "
        "eligibility is decided on pooled coherence at the selected support, so the "
        "profile-coherence values plotted here mark no eligibility decision, "
        "including for configured groups that were selected as supports. The lines "
        "define no validated quadrant. Profile coherence summarizes mean within-run "
        "sign homogeneity. It does not establish run-level sign recovery. "
        "Group-level uncertainty uses the zero-covariance approximation.\n"
    )


def build(paths: Paths) -> list[Path]:
    styles.apply_style()
    data, metadata_provenance = prepare_data(paths)
    serotypes = [s for s in SEROTYPE_ORDER if s in set(data["serotype"])]
    serotypes += sorted(set(data["serotype"].astype(str)) - set(serotypes))

    fig, (axa, axb) = plt.subplots(
        1, 2, figsize=(DOUBLE_COL, 3.7),
        gridspec_kw={"width_ratios": [1.25, 1.0], "wspace": 0.38})
    im, norm, cmap = _panel_heatmap(axa, data, serotypes)
    panel_letter(axa, "A")
    cbar = fig.colorbar(im, ax=axa, fraction=0.045, pad=0.03)
    # title on top of the colorbar instead of a side label, so it cannot
    # collide with Panel B's y-axis label
    cbar.ax.set_title("ρ", fontsize=styles.FS_TICK, pad=3)
    cbar.set_ticks([0, 0.5, 1.0])
    cbar.outline.set_linewidth(0.4)

    _panel_scatter(axb, data)
    panel_letter(axb, "B", dx=-0.08)

    handles, labels = axb.get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", bbox_to_anchor=(0.72, 0.005),
               fontsize=styles.FS_ANNOT, ncol=2, frameon=False)
    fig.text(0.72, 0.055, REFERENCE_NOTE,
             ha="center", va="bottom", fontsize=styles.FS_ANNOT, color="#7f7f7f")
    fig.suptitle(TITLE, x=0.5, y=0.995)
    # larger bottom margin keeps the x-axis label clear of the note and legend
    fig.subplots_adjust(bottom=0.24, top=0.88)
    written = save_figure(fig, paths, "figure3_catalytic_core")

    audit_rows = []
    for row in data.sort_values(["serotype", "chain", "domain"]).itertuples(index=False):
        audit_rows.append(
            {
                "dataset": row.serotype,
                "chain": row.chain,
                "configured_group": row.domain,
                "support_size": int(row.support_size),
                "rho": float(row.rho_domain),
                "coherence": float(row.coherence_domain),
                "is_configured_catalytic_group": bool(
                    row.is_configured_catalytic_group
                ),
            }
        )
    write_json(
        paths.output / "figure3_numeric_audit.json",
        {
            "configured_group_rows": len(data),
            "coherence_statistic": COHERENCE_STATISTIC,
            "metadata": metadata_provenance,
            "rows": audit_rows,
            "sources": {
                "outputs_s7/F3_domain_serotype_rho_heatmap.parquet": sha256(
                    paths.s7("F3_domain_serotype_rho_heatmap.parquet")
                ),
                "outputs_s7/F8_coherence_vs_rho.parquet": sha256(
                    paths.s7("F8_coherence_vs_rho.parquet")
                ),
            },
        },
    )
    (paths.output / "figure3_caption_draft.md").write_text(_caption_text())
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
