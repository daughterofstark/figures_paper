"""Supplementary Figure S2: conditional selected-support signed associations.

This presentation-only module consolidates the former main signed/BH plot and
the former sequence-axis supplement. It reads the frozen S4 selected-row screen
and the frozen S0 selected-support identifiers, then displays:

A. one fitted beta and implemented Wald interval per distinct finite fitted
   effect;
B. every selected item-row location, with copied fitted-effect representations
   explicitly outlined.

The stored row-level BH labels are displayed without recalculation. The figure
is conditional and descriptive: provisional support selection and coherence
eligibility precede the fitted estimates and implemented within-dataset BH
adjustment. No post-selection error-control claim is made.
"""
from __future__ import annotations

import hashlib
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.axes import Axes
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Rectangle

import styles
from figure_config import DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import panel_letter, parse_canon_label, save_figure

FIGURE_STEM = "Supplementary_Figure_S2"
DATASETS = tuple(SEROTYPE_ORDER)
AUTHORITATIVE_SOURCE_HASHES = {
    "outputs_s4/significance_screen.parquet": (
        "1a15af3264bfb9919cfd56d226bdb745450f53ecdbb02fe7967273f5d22a08e1"
    ),
    "outputs/stride_table.parquet": (
        "732509f98469cac7c5a5bb723ca78caa58341e47d4f5c1f2f6df3693fb89a213"
    ),
}
AUTHORITATIVE_COUNTS = {
    "selected_item_rows": 924,
    "distinct_selected_fits": 538,
    "finite_signed_rows": 302,
    "distinct_finite_fits": 268,
    "bh_flagged_rows": 152,
    "distinct_bh_flagged_fits": 145,
    "raw_interval_exclusion_rows": 181,
    "distinct_raw_interval_exclusion_fits": 174,
    "raw_wald_flagged_rows": 181,
}
AUTHORITATIVE_FAMILY_SIZES = {
    "DENV1": 70,
    "DENV2": 79,
    "DENV3": 94,
    "DENV4": 59,
}
EXPECTED_THRESHOLDS = {
    "provisional_rho_first_crossing": 0.5,
    "coherence_eligibility": 0.6,
}
FIT_KEY_COLUMNS = ("serotype", "mech_region_id")
STATUS_ORDER = (
    "ABSENT",
    "SELECTED_NO_FINITE",
    "NEGATIVE_UNFLAGGED",
    "POSITIVE_UNFLAGGED",
    "NEGATIVE_BH_FLAGGED",
    "POSITIVE_BH_FLAGGED",
)
STATUS_COLORS = {
    "ABSENT": "#FFFFFF",
    "SELECTED_NO_FINITE": "#D7DADD",
    "NEGATIVE_UNFLAGGED": "#B8CDD5",
    "POSITIVE_UNFLAGGED": "#D8BFC0",
    "NEGATIVE_BH_FLAGGED": "#355F72",
    "POSITIVE_BH_FLAGGED": "#8B5A5A",
}
FORBIDDEN_FIGURE_TERMS = (
    "statistically robust",
    "fdr-controlled",
    "significant mechanism",
    "discovery",
    "validated effect",
    "causal direction",
    "reproducible signed mechanism",
)

SCREEN_COLUMNS = {
    "serotype",
    "canon_label",
    "chain",
    "domain",
    "gated_scale_level",
    "direction",
    "is_signed",
    "beta_signed",
    "beta_se",
    "beta_ci_lower",
    "beta_ci_upper",
    "ci_excludes_zero",
    "p_value",
    "p_value_bh",
    "significant_raw",
    "significant_fdr",
}
STRIDE_COLUMNS = {
    "serotype",
    "canon_label",
    "scale_level",
    "is_gated_scale",
    "mech_region_id",
    "mech_n_loci",
    "mech_coherence",
    "gate_rho_star",
    "gate_coherence_threshold",
    "mech_beta_signed",
    "mech_beta_se",
}
FIT_EQUALITY_COLUMNS = (
    "scale_level",
    "mech_n_loci",
    "beta_signed",
    "beta_se",
    "beta_ci_lower",
    "beta_ci_upper",
    "ci_excludes_zero",
    "p_value",
    "p_value_bh",
    "significant_raw",
    "significant_fdr",
)


@dataclass(frozen=True)
class SignedFigureData:
    """Validated selected item rows and distinct fitted-effect records."""

    selected_rows: pd.DataFrame
    distinct_selected_fits: pd.DataFrame
    distinct_finite_fits: pd.DataFrame
    coordinate_order: pd.DataFrame
    source_hashes: Mapping[str, str]
    authoritative: bool


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_columns(frame: pd.DataFrame, required: set[str], source: Path) -> None:
    missing = sorted(required - set(frame.columns))
    if missing:
        raise ValueError(
            f"{source} is missing required columns {missing}; "
            f"available columns are {sorted(frame.columns)}"
        )


def _source_paths(paths: Paths) -> dict[str, Path]:
    return {
        "outputs_s4/significance_screen.parquet": paths.table(
            "s4", "significance_screen.parquet"
        ),
        "outputs/stride_table.parquet": paths.table("s0", "stride_table.parquet"),
    }


def _coordinate_order(rows: pd.DataFrame) -> pd.DataFrame:
    coordinates = rows[["canon_label", "chain"]].drop_duplicates().copy()
    if coordinates["canon_label"].duplicated().any():
        raise ValueError("a canonical coordinate maps to multiple chains")
    coordinates["_chain_order"] = coordinates["chain"].map(
        {chain: index for index, chain in enumerate(sorted(coordinates["chain"].unique()))}
    )
    coordinates["_residue"] = coordinates["canon_label"].map(
        lambda value: parse_canon_label(value)[1]
    )
    coordinates = coordinates.sort_values(
        ["_chain_order", "_residue", "canon_label"], kind="stable"
    ).reset_index(drop=True)
    coordinates["coordinate_index"] = np.arange(len(coordinates), dtype=int)
    return coordinates[
        ["coordinate_index", "canon_label", "chain", "_residue"]
    ].rename(columns={"_residue": "coordinate_number"})


def _fit_identity(rows: pd.DataFrame) -> pd.Series:
    return (
        rows["serotype"].astype(str)
        + "::"
        + rows["mech_region_id"].astype(str)
    )


def load_signed_figure_data(paths: Paths) -> SignedFigureData:
    """Load frozen selected rows and build an explicit fitted-effect mapping."""

    source_paths = _source_paths(paths)
    for path in source_paths.values():
        if not path.is_file():
            raise FileNotFoundError(f"missing required signed-figure source: {path}")
    source_hashes = {
        name: _sha256(path) for name, path in source_paths.items()
    }
    authoritative = source_hashes == AUTHORITATIVE_SOURCE_HASHES

    screen_path = source_paths["outputs_s4/significance_screen.parquet"]
    stride_path = source_paths["outputs/stride_table.parquet"]
    screen = pd.read_parquet(screen_path).copy()
    stride = pd.read_parquet(stride_path).copy()
    _require_columns(screen, SCREEN_COLUMNS, screen_path)
    _require_columns(stride, STRIDE_COLUMNS, stride_path)

    gated = stride[stride["is_gated_scale"].astype(bool)].copy()
    support_columns = [
        "serotype",
        "canon_label",
        "scale_level",
        "mech_region_id",
        "mech_n_loci",
        "mech_coherence",
        "gate_rho_star",
        "gate_coherence_threshold",
        "mech_beta_signed",
        "mech_beta_se",
    ]
    gated = gated[support_columns]
    if gated.duplicated(["serotype", "canon_label"]).any():
        raise ValueError("frozen selected-support table has multiple gated rows per item")
    if screen.duplicated(["serotype", "canon_label"]).any():
        raise ValueError("S4 screen has multiple selected rows per dataset-item")

    rows = screen.merge(
        gated,
        on=["serotype", "canon_label"],
        how="left",
        validate="one_to_one",
    )
    if rows["mech_region_id"].isna().any():
        missing = rows.loc[
            rows["mech_region_id"].isna(), ["serotype", "canon_label"]
        ]
        raise ValueError(
            "selected item rows lack explicit support identity: "
            f"{missing.head().to_dict(orient='records')}"
        )
    rows["fit_id"] = _fit_identity(rows)
    rows["finite_signed_test"] = np.isfinite(
        pd.to_numeric(rows["p_value"], errors="coerce")
    )
    rows["bh_flagged"] = rows["significant_fdr"].astype(bool)
    rows["raw_interval_exclusion"] = rows["ci_excludes_zero"].astype(bool)
    rows["raw_wald_flagged"] = rows["significant_raw"].astype(bool)

    multiplicity = rows.groupby("fit_id", sort=False).size().rename(
        "item_row_multiplicity"
    )
    rows = rows.join(multiplicity, on="fit_id")
    rows["repeated_representation"] = rows["item_row_multiplicity"] > 1

    finite_numeric_columns = (
        "beta_signed",
        "beta_se",
        "beta_ci_lower",
        "beta_ci_upper",
        "p_value",
        "p_value_bh",
    )
    for column in finite_numeric_columns:
        rows[column] = pd.to_numeric(rows[column], errors="coerce")

    for fit_id, group in rows.groupby("fit_id", sort=False):
        for column in FIT_EQUALITY_COLUMNS:
            if group[column].nunique(dropna=False) != 1:
                raise ValueError(
                    f"copied rows for fit {fit_id} disagree in {column}"
                )
    if (rows["bh_flagged"] & ~rows["finite_signed_test"]).any():
        raise ValueError("a BH-flagged row lacks a finite signed test")
    if not np.array_equal(
        rows["raw_interval_exclusion"].to_numpy(),
        rows["raw_wald_flagged"].to_numpy(),
    ):
        raise ValueError("raw interval exclusion and raw Wald flags disagree")
    if not np.array_equal(
        rows.loc[rows["finite_signed_test"], "bh_flagged"].to_numpy(),
        (
            rows.loc[rows["finite_signed_test"], "p_value_bh"].to_numpy(dtype=float)
            <= 0.05
        ),
    ):
        raise ValueError("stored BH labels do not match stored adjusted values")

    rows = rows.sort_values(
        ["serotype", "chain", "canon_label"], kind="stable"
    ).reset_index(drop=True)
    fits = (
        rows.sort_values(["serotype", "fit_id", "canon_label"], kind="stable")
        .drop_duplicates("fit_id", keep="first")
        .reset_index(drop=True)
    )
    finite_fits = (
        fits[fits["finite_signed_test"]]
        .sort_values(["serotype", "beta_signed", "fit_id"], kind="stable")
        .reset_index(drop=True)
    )
    coordinates = _coordinate_order(rows)

    data = SignedFigureData(
        selected_rows=rows,
        distinct_selected_fits=fits,
        distinct_finite_fits=finite_fits,
        coordinate_order=coordinates,
        source_hashes=source_hashes,
        authoritative=authoritative,
    )
    validate_signed_figure_data(data)
    return data


def bookkeeping_counts(data: SignedFigureData) -> dict[str, int]:
    """Return the locked row-level and distinct-fit bookkeeping counts."""

    rows = data.selected_rows
    fits = data.distinct_selected_fits
    return {
        "selected_item_rows": int(len(rows)),
        "distinct_selected_fits": int(len(fits)),
        "finite_signed_rows": int(rows["finite_signed_test"].sum()),
        "distinct_finite_fits": int(fits["finite_signed_test"].sum()),
        "bh_flagged_rows": int(rows["bh_flagged"].sum()),
        "distinct_bh_flagged_fits": int(fits["bh_flagged"].sum()),
        "raw_interval_exclusion_rows": int(
            rows["raw_interval_exclusion"].sum()
        ),
        "distinct_raw_interval_exclusion_fits": int(
            fits["raw_interval_exclusion"].sum()
        ),
        "raw_wald_flagged_rows": int(rows["raw_wald_flagged"].sum()),
    }


def validate_signed_figure_data(data: SignedFigureData) -> None:
    """Validate mapping, copied values, thresholds, and authoritative totals."""

    rows = data.selected_rows
    fits = data.distinct_selected_fits
    if rows["fit_id"].isna().any() or fits["fit_id"].duplicated().any():
        raise ValueError("fitted-effect identity is missing or non-unique")
    if not np.array_equal(
        rows["fit_id"].to_numpy(),
        _fit_identity(rows).to_numpy(),
    ):
        raise ValueError("fit_id does not equal dataset plus selected exact support")
    fit_counts = rows.groupby("fit_id", sort=False).size()
    observed_multiplicity = rows["fit_id"].map(fit_counts)
    if not np.array_equal(
        observed_multiplicity.to_numpy(dtype=int),
        rows["item_row_multiplicity"].to_numpy(dtype=int),
    ):
        raise ValueError("item-row multiplicity does not match explicit fit identity")

    for column in ("beta_signed", "beta_se", "beta_ci_lower", "beta_ci_upper"):
        finite = rows.loc[rows["finite_signed_test"], column]
        if not np.isfinite(finite.to_numpy(dtype=float)).all():
            raise ValueError(f"finite signed rows contain non-finite {column}")
    signs = np.sign(
        rows.loc[rows["finite_signed_test"], "beta_signed"].to_numpy(dtype=float)
    )
    if (signs == 0).any():
        raise ValueError("finite signed rows contain zero fitted beta")

    thresholds = rows[
        ["gate_rho_star", "gate_coherence_threshold"]
    ].drop_duplicates()
    if len(thresholds) != 1:
        raise ValueError("selected rows do not share one provisional rule pair")
    if not np.isclose(
        float(thresholds.iloc[0]["gate_rho_star"]),
        EXPECTED_THRESHOLDS["provisional_rho_first_crossing"],
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("selected rows do not use provisional rho reference 0.5")
    if not np.isclose(
        float(thresholds.iloc[0]["gate_coherence_threshold"]),
        EXPECTED_THRESHOLDS["coherence_eligibility"],
        rtol=0.0,
        atol=1e-12,
    ):
        raise ValueError("selected rows do not use coherence eligibility 0.6")

    if data.authoritative:
        counts = bookkeeping_counts(data)
        if counts != AUTHORITATIVE_COUNTS:
            raise ValueError(
                f"authoritative signed bookkeeping changed: {counts} "
                f"!= {AUTHORITATIVE_COUNTS}"
            )
        family_sizes = (
            rows.loc[rows["finite_signed_test"]]
            .groupby("serotype", sort=False)
            .size()
            .to_dict()
        )
        if family_sizes != AUTHORITATIVE_FAMILY_SIZES:
            raise ValueError(
                f"implemented within-dataset family sizes changed: {family_sizes}"
            )


def distinct_fits_by_dataset(
    data: SignedFigureData, dataset: str
) -> pd.DataFrame:
    """Return finite fitted effects in deterministic beta order."""

    return (
        data.distinct_finite_fits[
            data.distinct_finite_fits["serotype"] == dataset
        ]
        .sort_values(["beta_signed", "fit_id"], kind="stable")
        .reset_index(drop=True)
    )


def coordinate_status_matrix(
    data: SignedFigureData,
) -> tuple[np.ndarray, np.ndarray]:
    """Return row-level status codes and copied-representation mask."""

    coordinate_index = data.coordinate_order.set_index("canon_label")[
        "coordinate_index"
    ].to_dict()
    dataset_index = {dataset: index for index, dataset in enumerate(DATASETS)}
    matrix = np.zeros(
        (len(DATASETS), len(data.coordinate_order)), dtype=int
    )
    repeated = np.zeros_like(matrix, dtype=bool)
    status_codes = {status: index for index, status in enumerate(STATUS_ORDER)}

    for row in data.selected_rows.itertuples(index=False):
        if row.serotype not in dataset_index:
            continue
        i = dataset_index[str(row.serotype)]
        j = int(coordinate_index[str(row.canon_label)])
        if not bool(row.finite_signed_test):
            status = "SELECTED_NO_FINITE"
        else:
            sign = "POSITIVE" if float(row.beta_signed) > 0 else "NEGATIVE"
            suffix = "BH_FLAGGED" if bool(row.bh_flagged) else "UNFLAGGED"
            status = f"{sign}_{suffix}"
        matrix[i, j] = status_codes[status]
        repeated[i, j] = bool(
            row.repeated_representation and row.finite_signed_test
        )
    return matrix, repeated


def _panel_a(fig: Figure, spec, data: SignedFigureData) -> list[Axes]:
    inner = spec.subgridspec(1, 4, wspace=0.16)
    axes: list[Axes] = []
    y_limits = (
        float(data.distinct_finite_fits["beta_ci_lower"].min()) - 0.03,
        float(data.distinct_finite_fits["beta_ci_upper"].max()) + 0.03,
    )
    for index, dataset in enumerate(DATASETS):
        ax = fig.add_subplot(inner[0, index], sharey=axes[0] if axes else None)
        subset = distinct_fits_by_dataset(data, dataset)
        x = np.arange(len(subset))
        flagged = subset["bh_flagged"].to_numpy(dtype=bool)
        repeated = subset["repeated_representation"].to_numpy(dtype=bool)
        beta = subset["beta_signed"].to_numpy(dtype=float)
        lower = subset["beta_ci_lower"].to_numpy(dtype=float)
        upper = subset["beta_ci_upper"].to_numpy(dtype=float)

        ax.vlines(x, lower, upper, color="#B8B8B8", lw=0.45, alpha=0.72, zorder=1)
        ax.scatter(
            x[~flagged],
            beta[~flagged],
            s=10,
            color="#AEB3B7",
            edgecolor="white",
            linewidth=0.25,
            zorder=2,
            label="unflagged",
        )
        ax.scatter(
            x[flagged],
            beta[flagged],
            s=13,
            color="#2D3032",
            edgecolor="white",
            linewidth=0.25,
            zorder=3,
            label="BH-flagged",
        )
        if repeated.any():
            ax.scatter(
                x[repeated],
                beta[repeated],
                s=28,
                marker="s",
                facecolor="none",
                edgecolor=styles.OKABE_ITO["blue"],
                linewidth=0.65,
                zorder=4,
                label="represented by >1 item row",
            )
        ax.axhline(0, color="#555555", lw=0.55, zorder=0)
        ax.set_xlim(-1, max(len(subset), 1))
        ax.set_ylim(*y_limits)
        ax.set_xticks([])
        ax.set_title(
            f"{dataset}\n{len(subset)} fits / "
            f"{int(data.selected_rows.query('serotype == @dataset and finite_signed_test').shape[0])} rows",
            fontsize=styles.FS_LABEL,
            pad=5,
        )
        ax.set_xlabel("beta-ordered fits", fontsize=styles.FS_ANNOT)
        if index == 0:
            ax.set_ylabel("Fitted beta (implemented Wald interval)")
        else:
            ax.tick_params(axis="y", left=False, labelleft=False)
            ax.spines["left"].set_visible(False)
        axes.append(ax)
    return axes


def _panel_b(ax: Axes, data: SignedFigureData) -> None:
    matrix, repeated = coordinate_status_matrix(data)
    rows, columns = matrix.shape
    for i in range(rows):
        for j in range(columns):
            status = STATUS_ORDER[int(matrix[i, j])]
            rectangle = Rectangle(
                (j - 0.5, i - 0.5),
                1,
                1,
                facecolor=STATUS_COLORS[status],
                edgecolor="#E1E1E1",
                linewidth=0.12,
                hatch="////" if status.startswith("NEGATIVE") else None,
            )
            ax.add_patch(rectangle)
            if repeated[i, j]:
                ax.add_patch(
                    Rectangle(
                        (j - 0.5, i - 0.5),
                        1,
                        1,
                        facecolor="none",
                        edgecolor=styles.OKABE_ITO["blue"],
                        linewidth=0.42,
                    )
                )
    ax.set_xlim(-0.5, columns - 0.5)
    ax.set_ylim(rows - 0.5, -0.5)
    ax.set_yticks(range(rows), DATASETS)
    ax.tick_params(axis="y", length=0)

    coordinate_order = data.coordinate_order
    tick_indices = np.linspace(0, columns - 1, 6, dtype=int).tolist()
    tick_labels = (
        coordinate_order.set_index("coordinate_index")
        .loc[tick_indices, "canon_label"]
        .astype(str)
        .tolist()
    )
    ax.set_xticks(tick_indices, tick_labels, rotation=45, ha="right")
    for chain, subset in coordinate_order.groupby("chain", sort=False):
        start = int(subset["coordinate_index"].min())
        end = int(subset["coordinate_index"].max())
        ax.text(
            (start + end) / 2,
            1.04,
            str(chain),
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="bottom",
            fontsize=styles.FS_LABEL,
            fontweight="bold",
        )
        if start > 0:
            ax.axvline(start - 0.5, color="#222222", lw=0.7, zorder=5)

    ax.set_title(
        "Selected item-row locations and implemented BH status",
        loc="left",
        pad=13,
    )
    panel_letter(ax, "B", dx=-0.055, dy=1.12)
    handles = [
        Patch(
            facecolor=STATUS_COLORS["ABSENT"],
            edgecolor="#999999",
            label="Absent",
        ),
        Patch(
            facecolor=STATUS_COLORS["SELECTED_NO_FINITE"],
            edgecolor="none",
            label="Selected; no finite signed test",
        ),
        Patch(
            facecolor=STATUS_COLORS["NEGATIVE_UNFLAGGED"],
            edgecolor="#777777",
            hatch="////",
            label="Negative, unflagged",
        ),
        Patch(
            facecolor=STATUS_COLORS["POSITIVE_UNFLAGGED"],
            edgecolor="none",
            label="Positive, unflagged",
        ),
        Patch(
            facecolor=STATUS_COLORS["NEGATIVE_BH_FLAGGED"],
            edgecolor="#777777",
            hatch="////",
            label="Negative, BH-flagged",
        ),
        Patch(
            facecolor=STATUS_COLORS["POSITIVE_BH_FLAGGED"],
            edgecolor="none",
            label="Positive, BH-flagged",
        ),
        Patch(
            facecolor="none",
            edgecolor=styles.OKABE_ITO["blue"],
            linewidth=0.8,
            label="Copied fitted-effect representation",
        ),
    ]
    ax.legend(
        handles=handles,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.56),
        ncol=4,
        fontsize=styles.FS_ANNOT - 0.45,
        handlelength=1.1,
        columnspacing=0.8,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#AFAFAF")
        spine.set_linewidth(0.4)


def make_figure(data: SignedFigureData) -> Figure:
    """Construct the consolidated conditional signed-analysis supplement."""

    styles.apply_style()
    fig = plt.figure(figsize=(DOUBLE_COL, 6.1), constrained_layout=False)
    outer = fig.add_gridspec(
        2,
        1,
        height_ratios=(1.48, 0.64),
        left=0.085,
        right=0.985,
        bottom=0.25,
        top=0.84,
        hspace=0.60,
    )
    axes_a = _panel_a(fig, outer[0, 0], data)
    ax_b = fig.add_subplot(outer[1, 0])
    _panel_b(ax_b, data)
    panel_letter(axes_a[0], "A", dx=-0.20, dy=1.26)
    fig.text(
        0.085,
        0.885,
        "Distinct fitted beta estimates and implemented Wald intervals",
        ha="left",
        va="bottom",
        fontsize=styles.FS_TITLE,
    )
    fig.text(
        0.985,
        0.885,
        "Points: grey unflagged; black BH-flagged; blue square copied fit",
        ha="right",
        va="bottom",
        fontsize=styles.FS_ANNOT - 0.45,
        color="#555555",
    )
    fig.suptitle(
        "Conditional selected-support signed associations under the implemented provisional screen",
        x=0.085,
        y=0.975,
        ha="left",
        fontsize=styles.FS_TITLE,
        fontweight="bold",
    )
    return fig


def numeric_audit(data: SignedFigureData) -> dict[str, object]:
    """Return deterministic row/distinct-fit verification metadata."""

    rows = data.selected_rows
    fits = data.distinct_selected_fits
    finite_rows = rows[rows["finite_signed_test"]]
    counts = bookkeeping_counts(data)
    distinct_by_dataset = (
        data.distinct_finite_fits.groupby("serotype", sort=False)
        .size()
        .reindex(DATASETS, fill_value=0)
        .astype(int)
        .to_dict()
    )
    flagged_distinct_by_dataset = (
        fits[fits["bh_flagged"]]
        .groupby("serotype", sort=False)
        .size()
        .reindex(DATASETS, fill_value=0)
        .astype(int)
        .to_dict()
    )
    family_sizes = (
        finite_rows.groupby("serotype", sort=False)
        .size()
        .reindex(DATASETS, fill_value=0)
        .astype(int)
        .to_dict()
    )
    duplicate_fits = (
        fits[fits["repeated_representation"]][
            [
                "fit_id",
                "serotype",
                "scale_level",
                "mech_region_id",
                "mech_n_loci",
                "item_row_multiplicity",
                "finite_signed_test",
                "bh_flagged",
            ]
        ]
        .sort_values(["serotype", "fit_id"], kind="stable")
        .to_dict(orient="records")
    )
    return {
        "schema_version": "phase3i.supplementary_signed.audit.1",
        "interpretation": (
            "Conditional selected-support signed associations under the "
            "implemented provisional screen."
        ),
        "source_sha256": dict(data.source_hashes),
        "authoritative_frozen_sources": data.authoritative,
        "counts": counts,
        "within_dataset_finite_row_family_sizes": family_sizes,
        "distinct_finite_fits_by_dataset": distinct_by_dataset,
        "distinct_bh_flagged_fits_by_dataset": flagged_distinct_by_dataset,
        "distinct_signed_fits_represented_more_than_once": int(
            (
                fits["finite_signed_test"]
                & fits["repeated_representation"]
            ).sum()
        ),
        "repeated_finite_row_representations_beyond_first": int(
            (
                finite_rows.groupby("fit_id", sort=False).size() - 1
            ).clip(lower=0).sum()
        ),
        "duplicate_fit_records": duplicate_fits,
        "provisional_rho_reference": 0.5,
        "coherence_eligibility_reference": 0.6,
        "bh_family": "within_dataset_over_finite_selected_item_rows",
        "bh_values_recalculated": False,
        "post_selection_fdr_control_established": False,
        "wald_null_calibration_for_selected_pipeline_established": False,
    }


def _caption_text(data: SignedFigureData) -> str:
    counts = bookkeeping_counts(data)
    return (
        "# Supplementary Figure S2 caption draft\n\n"
        "**Supplementary Figure S2 | Conditional selected-support signed "
        "associations under the implemented provisional screen.** The frozen "
        f"screen contains {counts['selected_item_rows']} provisionally selected "
        f"item rows representing {counts['distinct_selected_fits']} distinct "
        "selected fitted effects. Signed testing follows data-dependent first-"
        "crossing selection at provisional rho=0.5 and data-dependent pooled "
        "coherence eligibility at 0.6. "
        f"(A) The {counts['finite_signed_rows']} finite signed-test rows represent "
        f"{counts['distinct_finite_fits']} distinct fitted effects; one fitted "
        "beta and its implemented Wald interval are shown per distinct fit. "
        "Square outlines mark fits copied to more than one selected item row. "
        f"The stored within-dataset BH procedure flagged "
        f"{counts['bh_flagged_rows']} item rows representing "
        f"{counts['distinct_bh_flagged_fits']} distinct fitted effects. "
        "(B) All selected item-row locations are shown. Darker cells carry the "
        "stored BH-adjusted flag, negative estimates are hatched, gray cells have "
        "no finite signed test, blank cells are absent, and blue outlines mark "
        "copied fitted-effect representations. BH adjustment was implemented "
        "separately within each dataset over finite selected item rows. Support "
        "selection and coherence eligibility reuse the fitted input values; "
        "post-selection false-discovery-rate control has not been established, "
        "and Wald null calibration has not been independently validated for this "
        "selected pipeline. The figure is conditional and descriptive.\n"
    )


def _write_metadata(paths: Paths, data: SignedFigureData) -> None:
    (paths.output / "Supplementary_Figure_S2_numeric_audit.json").write_text(
        json.dumps(numeric_audit(data), indent=2, sort_keys=True) + "\n"
    )
    (paths.output / "Supplementary_Figure_S2_source_provenance.json").write_text(
        json.dumps(
            {
                "schema_version": "phase3i.supplementary_signed.provenance.1",
                "sources": dict(data.source_hashes),
                "fit_identity": ["serotype", "mech_region_id"],
                "row_identity": ["serotype", "canon_label"],
                "plotting_module": Path(__file__).name,
                "stored_bh_labels_reused_without_recalculation": True,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    (paths.output / "Supplementary_Figure_S2_caption_draft.md").write_text(
        _caption_text(data)
    )


def build(paths: Paths) -> list[Path]:
    """Validate frozen selected rows, render S2, and write audit metadata."""

    data = load_signed_figure_data(paths)
    figure = make_figure(data)
    written = save_figure(figure, paths, FIGURE_STEM)
    _write_metadata(paths, data)
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for path in build(resolve_paths()):
        print(path)
