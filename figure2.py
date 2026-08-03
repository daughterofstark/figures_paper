"""Figure 2: support-specific magnitude stability changes with support size.

Scientific inputs are restricted to the deterministic Phase 3E presentation
tables in ``derived_data/phase3e``. The module verifies their hashes against the
Phase 3E manifest, validates the locked schemas and denominators, then renders:

A. observed item rho versus first unique configured-aggregate rho;
B. the corresponding observed delta-rho distributions;
C. completed V9A support-size summaries for separate zero-mean and fixed-mean
   zero-covariance synthetic characterizations;
D. observed first-crossing sensitivity at four provisional references, using
   collapsed unique supports only.

No raw/frozen table is read here and no STRIDE statistic is recomputed.
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

import styles
from figure_config import DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import panel_letter, save_figure

FIGURE_STEM = "figure2_spatial_resolution"
PHASE3E_MANIFEST_SHA256 = (
    "2043b2c418c5524964492ee32ea9cfad82393b2b8dd60db7b98f1690a6f89873"
)
EXPECTED_PROFILES = 924
EXPECTED_DATASETS = tuple(SEROTYPE_ORDER)
REGIME_ORDER = ("zero_mean", "frozen_scale_fixed_mean")
REGIME_LABELS = {
    "zero_mean": "Zero-mean benchmark",
    "frozen_scale_fixed_mean": "Fixed-mean benchmark",
}
SUPPORT_STRATA = ("1", "2-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129+")
CROSSING_CATEGORIES = (
    "item",
    "configured_group",
    "distinct_chain",
    "root",
    "unresolved",
)
CROSSING_LABELS = {
    "item": "Item",
    "configured_group": "Configured group",
    "distinct_chain": "Distinct chain",
    "root": "Whole-system root",
    "unresolved": "Unresolved",
}
CUTOFFS = (0.4, 0.5, 0.5330838393, 0.6)
FORBIDDEN_SCIENTIFIC_TERMS = (
    "rescue",
    "recovery",
    "emergence",
    "saturation",
    "supported scale",
    "licensed",
)

INPUT_FILES = (
    "figure2_item_first_aggregate.csv",
    "figure2_unique_supports.csv",
    "figure2_v9a_support_size.csv",
    "figure2_v9a_support_size_association.csv",
    "figure2_v9a_item_first_aggregate.csv",
    "figure2_v9a_first_crossing.csv",
    "figure2_v9a_alias_provenance.csv",
    "figure2_observed_first_crossing_sensitivity.csv",
)

REQUIRED_COLUMNS = {
    "figure2_item_first_aggregate.csv": {
        "dataset",
        "canonical_item_key",
        "item_rho",
        "first_aggregate_rho",
        "delta_rho",
        "first_aggregate_support_size",
    },
    "figure2_unique_supports.csv": {
        "dataset",
        "canonical_item_key",
        "collapsed_presentation_level",
        "support_size",
        "rho",
        "support_membership_hash",
    },
    "figure2_v9a_support_size.csv": {
        "dataset",
        "regime",
        "support_size_stratum",
        "support_size_stratum_index",
        "distinct_supports",
        "finite_fit_count",
        "mean_rho",
        "mc_t95_lower",
        "mc_t95_upper",
    },
    "figure2_v9a_support_size_association.csv": {
        "dataset",
        "regime",
        "mean_spearman",
        "mc_t95_lower",
        "mc_t95_upper",
        "median",
        "distribution_q025",
        "distribution_q975",
    },
    "figure2_v9a_item_first_aggregate.csv": {
        "dataset",
        "regime",
        "fraction_increased_mean",
        "fraction_decreased_mean",
        "median_delta_rho_mean",
    },
    "figure2_v9a_first_crossing.csv": {
        "dataset",
        "regime",
        "crossing_level",
        "crossing_level_index",
        "mean_draw_fraction",
    },
    "figure2_v9a_alias_provenance.csv": {"metric", "count", "fraction"},
    "figure2_observed_first_crossing_sensitivity.csv": {
        "cutoff",
        "dataset",
        "crossing_level",
        "crossing_level_index",
        "count",
        "denominator",
        "fraction",
    },
}


@dataclass(frozen=True)
class Figure2Data:
    """Validated Phase 3E tables used by Figure 2 and its metadata."""

    item_aggregate: pd.DataFrame
    unique_supports: pd.DataFrame
    v9a_support_size: pd.DataFrame
    v9a_association: pd.DataFrame
    v9a_item_aggregate: pd.DataFrame
    v9a_crossing: pd.DataFrame
    alias_provenance: pd.DataFrame
    observed_crossing: pd.DataFrame
    input_hashes: Mapping[str, str]


def default_data_dir() -> Path:
    """Return the repository-local Phase 3E presentation-data directory."""

    return Path(__file__).resolve().parent / "derived_data" / "phase3e"


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _require_columns(frame: pd.DataFrame, filename: str) -> None:
    missing = sorted(REQUIRED_COLUMNS[filename] - set(frame.columns))
    if missing:
        raise ValueError(
            f"{filename} is missing required columns {missing}; "
            f"available columns are {sorted(frame.columns)}"
        )


def load_phase3e_data(data_dir: Path | None = None) -> Figure2Data:
    """Load and validate only the eight locked Phase 3E Figure 2 tables."""

    root = Path(data_dir) if data_dir is not None else default_data_dir()
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"missing Phase 3E manifest: {manifest_path}")
    manifest_hash = _sha256(manifest_path)
    if manifest_hash != PHASE3E_MANIFEST_SHA256:
        raise ValueError(
            "Phase 3E manifest hash mismatch: "
            f"expected {PHASE3E_MANIFEST_SHA256}, observed {manifest_hash}"
        )
    manifest = json.loads(manifest_path.read_text())
    recorded = manifest.get("outputs", {})

    frames: dict[str, pd.DataFrame] = {}
    hashes: dict[str, str] = {}
    for filename in INPUT_FILES:
        path = root / filename
        if not path.is_file():
            raise FileNotFoundError(f"missing Phase 3E Figure 2 input: {path}")
        observed_hash = _sha256(path)
        expected_hash = recorded.get(filename, {}).get("sha256")
        if not expected_hash:
            raise ValueError(f"manifest has no SHA-256 for required input {filename}")
        if observed_hash != expected_hash:
            raise ValueError(
                f"hash mismatch for {filename}: expected {expected_hash}, "
                f"observed {observed_hash}"
            )
        frame = pd.read_csv(path)
        _require_columns(frame, filename)
        frames[filename] = frame
        hashes[filename] = observed_hash

    data = Figure2Data(
        item_aggregate=frames["figure2_item_first_aggregate.csv"],
        unique_supports=frames["figure2_unique_supports.csv"],
        v9a_support_size=frames["figure2_v9a_support_size.csv"],
        v9a_association=frames["figure2_v9a_support_size_association.csv"],
        v9a_item_aggregate=frames["figure2_v9a_item_first_aggregate.csv"],
        v9a_crossing=frames["figure2_v9a_first_crossing.csv"],
        alias_provenance=frames["figure2_v9a_alias_provenance.csv"],
        observed_crossing=frames[
            "figure2_observed_first_crossing_sensitivity.csv"
        ],
        input_hashes=hashes,
    )
    validate_figure2_data(data)
    return data


def validate_figure2_data(data: Figure2Data) -> None:
    """Fail loudly if locked Figure 2 denominators or categories drift."""

    paired = data.item_aggregate
    if len(paired) != EXPECTED_PROFILES:
        raise ValueError(f"expected {EXPECTED_PROFILES} profiles, found {len(paired)}")
    if paired.duplicated(["dataset", "canonical_item_key"]).any():
        raise ValueError("duplicate dataset-item rows in observed item/aggregate table")
    if set(paired["dataset"].astype(str)) != set(EXPECTED_DATASETS):
        raise ValueError("observed table does not contain exactly the four datasets")
    delta = paired["first_aggregate_rho"] - paired["item_rho"]
    if not np.allclose(delta, paired["delta_rho"], rtol=0.0, atol=1e-15):
        raise ValueError("stored delta_rho does not equal first_aggregate_rho-item_rho")

    unique = data.unique_supports
    unique_levels = set(unique["collapsed_presentation_level"].astype(str))
    expected_unique = {"item", "configured_group", "distinct_chain", "root"}
    if unique_levels != expected_unique:
        raise ValueError(
            f"unexpected collapsed support categories: {sorted(unique_levels)}"
        )
    if len(unique[unique["collapsed_presentation_level"] == "item"]) != EXPECTED_PROFILES:
        raise ValueError("unique-support table does not contain 924 singleton rows")

    support = data.v9a_support_size
    if set(support["regime"].astype(str)) != set(REGIME_ORDER):
        raise ValueError("V9A support-size regimes are not the locked two regimes")
    if set(support["dataset"].astype(str)) != set(EXPECTED_DATASETS):
        raise ValueError("V9A support-size table does not contain four datasets")
    if support.duplicated(["dataset", "regime", "support_size_stratum"]).any():
        raise ValueError("duplicate V9A dataset/regime/support-size rows")
    if not support["mc_t95_lower"].le(support["mean_rho"]).all():
        raise ValueError("V9A lower Monte Carlo interval exceeds mean rho")
    if not support["mean_rho"].le(support["mc_t95_upper"]).all():
        raise ValueError("V9A upper Monte Carlo interval is below mean rho")

    association = data.v9a_association
    if len(association) != 8 or association.duplicated(["dataset", "regime"]).any():
        raise ValueError("V9A association table must contain eight unique rows")

    v9a_aggregate = data.v9a_item_aggregate
    if len(v9a_aggregate) != 8 or v9a_aggregate.duplicated(
        ["dataset", "regime"]
    ).any():
        raise ValueError("V9A item/aggregate table must contain eight unique rows")

    v9a_crossing = data.v9a_crossing
    if set(v9a_crossing["crossing_level"].astype(str)) != set(CROSSING_CATEGORIES):
        raise ValueError("V9A crossing table contains non-collapsed categories")

    aliases = dict(
        zip(
            data.alias_provenance["metric"].astype(str),
            data.alias_provenance["count"].astype(int),
            strict=True,
        )
    )
    expected_aliases = {
        "total_decisions": 924000,
        "labels_changed": 672379,
        "selected_membership_changed": 0,
        "selected_rho_changed_gt_1e-12": 0,
    }
    if aliases != expected_aliases:
        raise ValueError(f"unexpected alias provenance: {aliases}")

    crossing = data.observed_crossing
    if set(crossing["crossing_level"].astype(str)) != set(CROSSING_CATEGORIES):
        raise ValueError("observed crossing table contains non-collapsed categories")
    observed_cutoffs = sorted(crossing["cutoff"].astype(float).unique())
    if len(observed_cutoffs) != 4 or not np.allclose(
        observed_cutoffs, CUTOFFS, rtol=0.0, atol=1e-12
    ):
        raise ValueError(f"unexpected provisional references: {observed_cutoffs}")
    overall = crossing[crossing["dataset"] == "ALL"]
    if len(overall) != len(CUTOFFS) * len(CROSSING_CATEGORIES):
        raise ValueError("overall crossing sensitivity must have 20 rows")
    if not overall["denominator"].eq(EXPECTED_PROFILES).all():
        raise ValueError("overall crossing-sensitivity denominator is not 924")
    totals = overall.groupby("cutoff", sort=True)["count"].sum()
    if not totals.eq(EXPECTED_PROFILES).all():
        raise ValueError("crossing categories do not sum to 924 at every cutoff")


def _dataset_order(frame: pd.DataFrame) -> list[str]:
    present = set(frame["dataset"].astype(str))
    return [dataset for dataset in EXPECTED_DATASETS if dataset in present]


def _panel_a(ax: Axes, data: Figure2Data) -> None:
    paired = data.item_aggregate
    for dataset in _dataset_order(paired):
        subset = paired[paired["dataset"] == dataset]
        ax.scatter(
            subset["item_rho"],
            subset["first_aggregate_rho"],
            s=8,
            alpha=0.42,
            color=styles.serotype_color(dataset),
            edgecolors="none",
            rasterized=False,
            label=f"{dataset} (n={len(subset)})",
            zorder=2,
        )
    ax.plot([0, 1], [0, 1], color="#5F6368", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set(xlim=(-0.02, 1.02), ylim=(-0.02, 1.02))
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel(r"Item $\rho$")
    ax.set_ylabel(r"First unique aggregate $\rho$")
    ax.set_title("Observed item and first aggregate", loc="left")
    ax.text(
        0.03,
        0.04,
        f"n = {len(paired)}",
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=styles.FS_ANNOT,
        color="#444444",
    )
    ax.legend(
        loc="lower right",
        ncol=2,
        fontsize=styles.FS_ANNOT - 0.3,
        handletextpad=0.35,
        columnspacing=0.7,
    )


def _panel_b(ax: Axes, data: Figure2Data) -> None:
    paired = data.item_aggregate
    datasets = _dataset_order(paired)
    values = [
        paired.loc[paired["dataset"] == dataset, "delta_rho"].to_numpy(dtype=float)
        for dataset in datasets
    ]
    parts = ax.violinplot(
        values,
        positions=np.arange(len(datasets)),
        widths=0.78,
        showmeans=False,
        showmedians=False,
        showextrema=False,
        bw_method=0.25,
    )
    for body, dataset in zip(parts["bodies"], datasets, strict=True):
        body.set_facecolor(styles.serotype_color(dataset))
        body.set_edgecolor(styles.serotype_color(dataset))
        body.set_alpha(0.22)
        body.set_linewidth(0.6)

    for index, (dataset, vals) in enumerate(zip(datasets, values, strict=True)):
        q1, median, q3 = np.quantile(vals, [0.25, 0.5, 0.75])
        ax.vlines(index, q1, q3, color=styles.serotype_color(dataset), lw=4.0, zorder=3)
        ax.plot(
            index,
            median,
            marker="_",
            markersize=8,
            markeredgewidth=1.2,
            color="#111111",
            zorder=4,
        )
        positive = np.mean(vals > 0)
        negative = np.mean(vals < 0)
        ax.text(
            index,
            0.985,
            f"n={len(vals)}\n+ {positive:.1%}\n- {negative:.1%}",
            transform=ax.get_xaxis_transform(),
            ha="center",
            va="top",
            fontsize=styles.FS_ANNOT - 0.4,
            color="#3F3F3F",
            linespacing=1.05,
        )

    ax.axhline(0, color="#3F3F3F", lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.set_xticks(range(len(datasets)), datasets)
    ax.set_ylabel(r"$\Delta\rho$ (first aggregate - item)")
    ax.set_title("Observed change by dataset", loc="left")


def _plot_v9a_regime(ax: Axes, data: Figure2Data, regime: str) -> None:
    support = data.v9a_support_size
    association = data.v9a_association
    subset = support[support["regime"] == regime]
    for dataset in EXPECTED_DATASETS:
        rows = subset[subset["dataset"] == dataset].sort_values(
            "support_size_stratum_index", kind="stable"
        )
        x = rows["support_size_stratum_index"].to_numpy(dtype=float)
        y = rows["mean_rho"].to_numpy(dtype=float)
        lower = rows["mc_t95_lower"].to_numpy(dtype=float)
        upper = rows["mc_t95_upper"].to_numpy(dtype=float)
        ax.errorbar(
            x,
            y,
            yerr=np.vstack((y - lower, upper - y)),
            color=styles.serotype_color(dataset),
            marker="o",
            markersize=2.5,
            lw=0.8,
            elinewidth=0.6,
            capsize=1.3,
            label=dataset,
            zorder=2,
        )
    ax.set_xlim(-0.3, len(SUPPORT_STRATA) - 0.7)
    ax.set_ylim(0, 1.02)
    ax.set_xticks(range(len(SUPPORT_STRATA)))
    ax.set_xticklabels(SUPPORT_STRATA, rotation=45, ha="right")
    ax.set_title(REGIME_LABELS[regime], fontsize=styles.FS_LABEL, loc="left", pad=2)
    ax.set_xlabel("Support size")

    assoc = association[association["regime"] == regime].set_index("dataset")
    values = [float(assoc.loc[dataset, "mean_spearman"]) for dataset in EXPECTED_DATASETS]
    annotation = (
        "Mean Spearman\n"
        + "  ".join(f"D{i + 1} {value:.3f}" for i, value in enumerate(values[:2]))
        + "\n"
        + "  ".join(f"D{i + 3} {value:.3f}" for i, value in enumerate(values[2:]))
    )
    ax.text(
        0.03,
        0.04,
        annotation,
        transform=ax.transAxes,
        ha="left",
        va="bottom",
        fontsize=styles.FS_ANNOT - 0.6,
        color="#333333",
        linespacing=1.08,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 1.2},
    )


def _panel_d(ax: Axes, data: Figure2Data) -> None:
    overall = data.observed_crossing[
        data.observed_crossing["dataset"] == "ALL"
    ].copy()
    overall["crossing_level"] = pd.Categorical(
        overall["crossing_level"], categories=CROSSING_CATEGORIES, ordered=True
    )
    overall = overall.sort_values(["cutoff", "crossing_level"], kind="stable")

    colors = {
        "item": "#3B6C8E",
        "configured_group": "#79A9A4",
        "distinct_chain": "#B8B8B8",
        "root": "#D7D7D7",
        "unresolved": "#F1F1F1",
    }
    y = np.arange(len(CUTOFFS))
    left = np.zeros(len(CUTOFFS), dtype=float)
    for category in CROSSING_CATEGORIES:
        rows = overall[overall["crossing_level"] == category].sort_values("cutoff")
        fractions = rows["fraction"].to_numpy(dtype=float)
        counts = rows["count"].to_numpy(dtype=int)
        ax.barh(
            y,
            fractions,
            left=left,
            height=0.62,
            color=colors[category],
            edgecolor="white",
            linewidth=0.45,
            label=CROSSING_LABELS[category],
        )
        for row_index, (start, fraction, count) in enumerate(
            zip(left, fractions, counts, strict=True)
        ):
            if fraction >= 0.065:
                text_color = "white" if category == "item" else "#222222"
                ax.text(
                    start + fraction / 2,
                    row_index,
                    str(count),
                    ha="center",
                    va="center",
                    fontsize=styles.FS_ANNOT - 0.4,
                    color=text_color,
                )
        left += fractions

    ax.set_xlim(0, 1)
    ax.set_xticks([0, 0.5, 1], ["0", "50", "100"])
    ax.set_yticks(y, ["0.400", "0.500", "0.533", "0.600"])
    ax.invert_yaxis()
    ax.set_xlabel("Profiles (%)")
    ax.set_ylabel(r"Provisional $\rho$ reference")
    ax.set_title("Provisional first crossing (unique supports)", loc="left")
    ax.legend(
        loc="upper center",
        bbox_to_anchor=(0.5, -0.23),
        ncol=3,
        fontsize=styles.FS_ANNOT - 0.5,
        handlelength=1.0,
        columnspacing=0.8,
    )


def make_figure(data: Figure2Data) -> Figure:
    """Construct the locked four-panel Figure 2 from validated data."""

    styles.apply_style()
    fig = plt.figure(figsize=(DOUBLE_COL, 6.05), constrained_layout=False)
    outer = fig.add_gridspec(
        2,
        2,
        width_ratios=(1.0, 1.14),
        height_ratios=(1.0, 0.92),
        left=0.075,
        right=0.985,
        bottom=0.115,
        top=0.93,
        wspace=0.27,
        hspace=0.42,
    )

    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    _panel_a(ax_a, data)
    _panel_b(ax_b, data)
    panel_letter(ax_a, "A", dx=-0.08, dy=1.02)
    panel_letter(ax_b, "B", dx=-0.07, dy=1.02)

    panel_c = outer[1, 0].subgridspec(1, 2, wspace=0.10)
    ax_c0 = fig.add_subplot(panel_c[0, 0])
    ax_c1 = fig.add_subplot(panel_c[0, 1], sharey=ax_c0)
    _plot_v9a_regime(ax_c0, data, "zero_mean")
    _plot_v9a_regime(ax_c1, data, "frozen_scale_fixed_mean")
    ax_c0.set_ylabel(r"Mean $\rho$ across draws")
    ax_c1.tick_params(labelleft=False)
    ax_c1.spines["left"].set_visible(False)
    ax_c0.text(
        0,
        1.13,
        "V9A: zero-covariance synthetic characterization",
        transform=ax_c0.transAxes,
        fontsize=styles.FS_TITLE,
        fontweight="normal",
        ha="left",
        va="bottom",
    )
    panel_letter(ax_c0, "C", dx=-0.16, dy=1.10)

    ax_d = fig.add_subplot(outer[1, 1])
    _panel_d(ax_d, data)
    panel_letter(ax_d, "D", dx=-0.07, dy=1.02)

    fig.suptitle(
        "Support-specific magnitude stability changes with support size",
        x=0.075,
        y=0.982,
        ha="left",
        fontsize=styles.FS_TITLE,
        fontweight="bold",
    )
    return fig


def numeric_audit(data: Figure2Data) -> dict[str, object]:
    """Return deterministic panel-verification metadata."""

    paired = data.item_aggregate
    panel_b: list[dict[str, object]] = []
    for dataset in _dataset_order(paired):
        subset = paired[paired["dataset"] == dataset]
        delta = subset["delta_rho"].to_numpy(dtype=float)
        panel_b.append(
            {
                "dataset": dataset,
                "n": int(len(subset)),
                "median_delta_rho": float(np.median(delta)),
                "fraction_positive": float(np.mean(delta > 0)),
                "fraction_negative": float(np.mean(delta < 0)),
            }
        )

    association = data.v9a_association.sort_values(
        ["regime", "dataset"], kind="stable"
    )
    panel_c_association = [
        {
            "dataset": str(row.dataset),
            "regime": str(row.regime),
            "mean_spearman": float(row.mean_spearman),
            "mc_t95_lower": float(row.mc_t95_lower),
            "mc_t95_upper": float(row.mc_t95_upper),
            "median": float(row.median),
            "distribution_q025": float(row.distribution_q025),
            "distribution_q975": float(row.distribution_q975),
        }
        for row in association.itertuples(index=False)
    ]
    support_strata = {}
    for regime in REGIME_ORDER:
        present = set(
            data.v9a_support_size.loc[
                data.v9a_support_size["regime"] == regime,
                "support_size_stratum",
            ].astype(str)
        )
        support_strata[regime] = [
            stratum for stratum in SUPPORT_STRATA if stratum in present
        ]

    crossing = data.observed_crossing[
        data.observed_crossing["dataset"] == "ALL"
    ].sort_values(["cutoff", "crossing_level_index"], kind="stable")
    panel_d = [
        {
            "cutoff": float(row.cutoff),
            "category": str(row.crossing_level),
            "count": int(row.count),
            "denominator": int(row.denominator),
            "fraction": float(row.fraction),
        }
        for row in crossing.itertuples(index=False)
    ]

    aliases = {
        str(row.metric): {
            "count": int(row.count),
            "fraction": float(row.fraction),
        }
        for row in data.alias_provenance.itertuples(index=False)
    }
    zero_mean_aggregate = data.v9a_item_aggregate[
        data.v9a_item_aggregate["regime"] == "zero_mean"
    ].sort_values("dataset", kind="stable")

    return {
        "schema_version": "phase3f.figure2.audit.1",
        "scientific_inputs": list(INPUT_FILES),
        "input_sha256": dict(sorted(data.input_hashes.items())),
        "panel_a": {
            "n": int(len(paired)),
            "fraction_y_gt_x": float(
                np.mean(paired["first_aggregate_rho"] > paired["item_rho"])
            ),
            "fraction_y_lt_x": float(
                np.mean(paired["first_aggregate_rho"] < paired["item_rho"])
            ),
            "fraction_y_eq_x": float(
                np.mean(paired["first_aggregate_rho"] == paired["item_rho"])
            ),
            "median_item_rho": float(paired["item_rho"].median()),
            "median_first_aggregate_rho": float(
                paired["first_aggregate_rho"].median()
            ),
        },
        "panel_b": panel_b,
        "panel_c": {
            "support_size_strata": support_strata,
            "association": panel_c_association,
            "zero_mean_fraction_increased": [
                {
                    "dataset": str(row.dataset),
                    "fraction": float(row.fraction_increased_mean),
                }
                for row in zero_mean_aggregate.itertuples(index=False)
            ],
        },
        "panel_d": panel_d,
        "alias_provenance": aliases,
    }


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_metadata(data: Figure2Data, output_dir: Path) -> list[Path]:
    """Write deterministic Figure 2 verification, provenance, and caption files."""

    audit = numeric_audit(data)
    audit_path = output_dir / "figure2_numeric_audit.json"
    _write_json(audit_path, audit)

    provenance_path = output_dir / "figure2_source_provenance.json"
    _write_json(
        provenance_path,
        {
            "schema_version": "phase3f.figure2.provenance.1",
            "phase3e_manifest_sha256": PHASE3E_MANIFEST_SHA256,
            "scientific_inputs": [
                {"path": filename, "sha256": data.input_hashes[filename]}
                for filename in INPUT_FILES
            ],
        },
    )

    aliases = audit["alias_provenance"]
    zero_fractions = [
        row["fraction"] for row in audit["panel_c"]["zero_mean_fraction_increased"]
    ]
    caption = (
        "# Figure 2 caption draft\n\n"
        "**Figure 2 | Support-specific magnitude stability changes with support "
        "size.** (A) Observed item-level magnitude-stability statistics are "
        "compared with the corresponding first distinct configured-aggregate "
        f"statistics for {EXPECTED_PROFILES} frozen dataset-item profiles; the "
        "diagonal marks equality. (B) Observed first-aggregate-minus-item "
        "differences are shown descriptively by dataset. Exact support aliases "
        "were collapsed by membership before presentation. (C) Completed V9A "
        "K=3 zero-covariance synthetic characterizations show mean rho by support-"
        "size stratum, with 95% Monte Carlo intervals for draw-level means. "
        "Zero-mean item-to-first-aggregate increases exceeded 99.99% in every "
        f"dataset template (minimum {min(zero_fractions):.4%}), and mean draw-level "
        "Spearman coefficients were positive under zero mean. V9A is a synthetic "
        "characterization, not a model of the upstream process. These results show "
        "that an aggregate increase is not evidence of recovered signal and that "
        "rho values at different supports describe support-specific objects. "
        "(D) First crossings across unique supports are shown only as sensitivity "
        "to four provisional operational references. V9A conditions on alias "
        "equivalence classes verified separately from frozen production "
        "memberships. Collapsed traversal is exact-support deduplication of named "
        "traversal, so equal selection follows structurally once duplicate "
        "representations share rho; it is not a Monte Carlo invariance rate. "
        "Production rho was equal across all 72,000 duplicate-representation "
        f"checks. Collapsing aliases changed {aliases['labels_changed']['count']:,} "
        f"of {aliases['total_decisions']['count']:,} reported item-decision labels.\n"
    )
    caption_path = output_dir / "figure2_caption_draft.md"
    caption_path.write_text(caption)
    return [audit_path, provenance_path, caption_path]


def build(
    paths: Paths,
    *,
    data_dir: Path | None = None,
    write_supporting_metadata: bool = True,
) -> list[Path]:
    """Render Figure 2 from verified Phase 3E inputs."""

    data = load_phase3e_data(data_dir)
    fig = make_figure(data)
    written = save_figure(fig, paths, FIGURE_STEM)
    if write_supporting_metadata:
        write_metadata(data, paths.output)
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for output_path in build(resolve_paths()):
        print(output_path)
