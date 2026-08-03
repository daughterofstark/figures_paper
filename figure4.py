"""Figure 4: single-run top-item rankings vary across replicate runs.

Scientific inputs are restricted to the deterministic Phase 3E presentation
tables in ``derived_data/phase3e``. The module verifies their hashes against the
Phase 3E manifest, validates the locked rankings and overlap records, then shows:

A. the item ranked first by absolute effect in each frozen run;
B. all within-dataset pairwise top-N Jaccard overlaps;
C. rank, membership, and sign for the union of top-10 items across runs.

No raw replicate table is read here. Ranking and Jaccard values are never
recomputed by the plotting layer.
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
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap
from matplotlib.figure import Figure
from matplotlib.lines import Line2D

import styles
from figure_config import DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import panel_letter, save_figure

FIGURE_STEM = "figure4_single_run_rankings"
PHASE3E_MANIFEST_SHA256 = (
    "2043b2c418c5524964492ee32ea9cfad82393b2b8dd60db7b98f1690a6f89873"
)
INPUT_FILES = (
    "figure4_top1.csv",
    "figure4_topn_jaccard.csv",
    "figure4_top10_heatmap.csv",
)
DATASETS = tuple(SEROTYPE_ORDER)
RUN_INDICES = (1, 2, 3)
N_VALUES = (5, 10, 20)
EXPECTED_TOP1_ROWS = 12
EXPECTED_JACCARD_ROWS = 36
EXPECTED_HEATMAP_ROWS = 315
SIGN_COLORS = {1: "#355F72", -1: "#8B5A5A", 0: "#777777"}
SIGN_MARKERS = {1: "^", -1: "v", 0: "o"}
PAIR_MARKERS = {(1, 2): "o", (1, 3): "s", (2, 3): "^"}
FORBIDDEN_FIGURE_TERMS = (
    "supportable",
    "supported scale",
    "recovery",
    "recovered",
    "mechanism",
    "significant",
    "robust",
    "true driver",
)

REQUIRED_COLUMNS = {
    "figure4_top1.csv": {
        "dataset",
        "run",
        "run_index",
        "top_item",
        "canonical_coordinate",
        "effect",
        "sign",
        "absolute_effect",
    },
    "figure4_topn_jaccard.csv": {
        "dataset",
        "N",
        "run_a",
        "run_a_index",
        "run_b",
        "run_b_index",
        "intersection_count",
        "union_count",
        "jaccard",
    },
    "figure4_top10_heatmap.csv": {
        "dataset",
        "union_order",
        "canonical_coordinate",
        "run",
        "run_index",
        "rank",
        "effect",
        "sign",
        "absolute_effect",
        "in_top10",
    },
}


@dataclass(frozen=True)
class Figure4Data:
    """Validated Phase 3E tables used by Figure 4."""

    top1: pd.DataFrame
    jaccard: pd.DataFrame
    heatmap: pd.DataFrame
    input_hashes: Mapping[str, str]
    tie_count: int
    tied_rows: int
    tie_rule: str


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


def load_phase3e_data(data_dir: Path | None = None) -> Figure4Data:
    """Load and validate only the three locked Phase 3E Figure 4 tables."""

    root = Path(data_dir) if data_dir is not None else default_data_dir()
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"missing Phase 3E manifest: {manifest_path}")
    observed_manifest_hash = _sha256(manifest_path)
    if observed_manifest_hash != PHASE3E_MANIFEST_SHA256:
        raise ValueError(
            "Phase 3E manifest hash mismatch: "
            f"expected {PHASE3E_MANIFEST_SHA256}, observed {observed_manifest_hash}"
        )
    manifest = json.loads(manifest_path.read_text())
    recorded = manifest.get("outputs", {})

    frames: dict[str, pd.DataFrame] = {}
    hashes: dict[str, str] = {}
    for filename in INPUT_FILES:
        path = root / filename
        if not path.is_file():
            raise FileNotFoundError(f"missing Phase 3E Figure 4 input: {path}")
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

    ranking = manifest.get("ranking", {})
    if ranking.get("primary_order") != "absolute_effect_descending":
        raise ValueError("Phase 3E ranking provenance is not absolute-effect order")
    data = Figure4Data(
        top1=frames["figure4_top1.csv"],
        jaccard=frames["figure4_topn_jaccard.csv"],
        heatmap=frames["figure4_top10_heatmap.csv"],
        input_hashes=hashes,
        tie_count=int(ranking.get("exact_tie_groups", -1)),
        tied_rows=int(ranking.get("rows_in_exact_ties", -1)),
        tie_rule=str(ranking.get("tie_break", "")),
    )
    validate_figure4_data(data)
    return data


def validate_figure4_data(data: Figure4Data) -> None:
    """Fail loudly if locked ranking records, ordering, or arithmetic drift."""

    top1 = data.top1
    if len(top1) != EXPECTED_TOP1_ROWS:
        raise ValueError(f"expected 12 top-1 rows, found {len(top1)}")
    if top1.duplicated(["dataset", "run_index"]).any():
        raise ValueError("duplicate dataset-run row in top-1 table")
    expected_runs = pd.MultiIndex.from_product(
        [DATASETS, RUN_INDICES], names=["dataset", "run_index"]
    )
    observed_runs = pd.MultiIndex.from_frame(top1[["dataset", "run_index"]])
    if set(observed_runs) != set(expected_runs):
        raise ValueError("top-1 table does not contain four datasets x three runs")
    if not np.allclose(
        top1["absolute_effect"],
        np.abs(top1["effect"]),
        rtol=0.0,
        atol=1e-15,
    ):
        raise ValueError("top-1 absolute effect does not equal abs(effect)")
    if not np.array_equal(
        top1["sign"].to_numpy(dtype=int),
        np.sign(top1["effect"]).astype(int).to_numpy(),
    ):
        raise ValueError("top-1 sign does not match retained effect sign")

    jaccard = data.jaccard
    if len(jaccard) != EXPECTED_JACCARD_ROWS:
        raise ValueError(f"expected 36 Jaccard rows, found {len(jaccard)}")
    if set(jaccard["dataset"].astype(str)) != set(DATASETS):
        raise ValueError("Jaccard table does not contain exactly four datasets")
    if set(jaccard["N"].astype(int)) != set(N_VALUES):
        raise ValueError("Jaccard table does not contain exactly N=5,10,20")
    pair_columns = ["dataset", "N", "run_a_index", "run_b_index"]
    if jaccard.duplicated(pair_columns).any():
        raise ValueError("duplicate dataset-N-run-pair Jaccard row")
    expected_pairs = set(PAIR_MARKERS)
    for (dataset, n_value), subset in jaccard.groupby(["dataset", "N"], sort=False):
        pairs = set(
            zip(
                subset["run_a_index"].astype(int),
                subset["run_b_index"].astype(int),
                strict=True,
            )
        )
        if pairs != expected_pairs:
            raise ValueError(
                f"{dataset}, N={n_value} does not contain exactly three run pairs"
            )
    if not jaccard["jaccard"].between(0.0, 1.0).all():
        raise ValueError("Jaccard values lie outside [0,1]")
    quotient = jaccard["intersection_count"] / jaccard["union_count"]
    if not np.allclose(quotient, jaccard["jaccard"], rtol=0.0, atol=1e-15):
        raise ValueError("stored Jaccard does not match intersection/union")

    heatmap = data.heatmap
    if len(heatmap) != EXPECTED_HEATMAP_ROWS:
        raise ValueError(f"expected 315 heatmap rows, found {len(heatmap)}")
    if heatmap.duplicated(["dataset", "union_order", "run_index"]).any():
        raise ValueError("duplicate dataset-union-order-run heatmap row")
    for dataset in DATASETS:
        subset = heatmap[heatmap["dataset"] == dataset]
        if subset.empty:
            raise ValueError(f"heatmap is missing dataset {dataset}")
        if list(subset[["union_order", "run_index"]].itertuples(index=False, name=None)) != sorted(
            subset[["union_order", "run_index"]].itertuples(index=False, name=None)
        ):
            raise ValueError(f"{dataset} heatmap rows do not retain Phase 3E order")
        per_item = subset.groupby(
            ["union_order", "canonical_coordinate"], sort=False
        )["run_index"].agg(list)
        if not all(runs == [1, 2, 3] for runs in per_item):
            raise ValueError(f"{dataset} heatmap item lacks the three ordered runs")
        expected_orders = list(range(int(subset["union_order"].max()) + 1))
        if list(per_item.index.get_level_values("union_order")) != expected_orders:
            raise ValueError(f"{dataset} union_order is not contiguous and deterministic")
    if not (heatmap["rank"].astype(int) >= 1).all():
        raise ValueError("heatmap contains a non-positive rank")
    expected_membership = heatmap["rank"].astype(int) <= 10
    if not np.array_equal(expected_membership, heatmap["in_top10"].astype(bool)):
        raise ValueError("explicit top-10 membership does not match retained rank")
    if not np.array_equal(
        heatmap["sign"].to_numpy(dtype=int),
        np.sign(heatmap["effect"]).astype(int).to_numpy(),
    ):
        raise ValueError("heatmap sign does not match retained effect sign")

    forbidden_columns = {"rho", "cutoff", "threshold", "p_value", "significant"}
    for name, frame in (
        ("top1", top1),
        ("jaccard", jaccard),
        ("heatmap", heatmap),
    ):
        leaked = forbidden_columns.intersection(frame.columns)
        if leaked:
            raise ValueError(f"{name} table contains forbidden columns {sorted(leaked)}")
    if data.tie_count < 0 or data.tied_rows < 0 or not data.tie_rule:
        raise ValueError("ranking tie provenance is incomplete")


def top1_unique_counts(data: Figure4Data) -> dict[str, int]:
    """Return distinct top-1 nominations per dataset."""

    return {
        dataset: int(
            data.top1.loc[
                data.top1["dataset"] == dataset, "canonical_coordinate"
            ].nunique()
        )
        for dataset in DATASETS
    }


def jaccard_ranges(data: Figure4Data) -> dict[str, dict[str, float]]:
    """Return descriptive min/median/max values by dataset and N."""

    ranges: dict[str, dict[str, float]] = {}
    for (dataset, n_value), subset in data.jaccard.groupby(
        ["dataset", "N"], sort=False
    ):
        values = subset["jaccard"].to_numpy(dtype=float)
        ranges[f"{dataset}|{int(n_value)}"] = {
            "min": float(np.min(values)),
            "median": float(np.median(values)),
            "max": float(np.max(values)),
        }
    return ranges


def heatmap_matrix(
    data: Figure4Data, dataset: str
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Return top-10 ranks, signs, and item labels in supplied union order."""

    subset = data.heatmap[data.heatmap["dataset"] == dataset]
    item_rows = (
        subset[["union_order", "canonical_coordinate"]]
        .drop_duplicates()
        .sort_values("union_order", kind="stable")
    )
    labels = item_rows["canonical_coordinate"].astype(str).tolist()
    ranks = np.full((len(labels), 3), np.nan)
    signs = np.zeros((len(labels), 3), dtype=int)
    for row in subset.itertuples(index=False):
        i = int(row.union_order)
        j = int(row.run_index) - 1
        if bool(row.in_top10):
            ranks[i, j] = int(row.rank)
            signs[i, j] = int(row.sign)
    return ranks, signs, labels


def _panel_a(ax: Axes, data: Figure4Data) -> None:
    top1 = data.top1.set_index(["dataset", "run_index"])
    for y_index, dataset in enumerate(DATASETS):
        for run_index in RUN_INDICES:
            row = top1.loc[(dataset, run_index)]
            sign = int(row["sign"])
            ax.scatter(
                run_index,
                y_index - 0.12,
                marker=SIGN_MARKERS[sign],
                s=24,
                facecolor=SIGN_COLORS[sign],
                edgecolor="white",
                linewidth=0.35,
                zorder=3,
            )
            ax.text(
                run_index,
                y_index + 0.14,
                str(row["canonical_coordinate"]),
                ha="center",
                va="center",
                fontsize=styles.FS_ANNOT,
                color="#222222",
            )
    for x in (1.5, 2.5):
        ax.axvline(x, color="#E2E2E2", lw=0.45, zorder=0)
    for y in (0.5, 1.5, 2.5):
        ax.axhline(y, color="#E2E2E2", lw=0.45, zorder=0)
    ax.set_xlim(0.5, 3.5)
    ax.set_ylim(3.55, -0.55)
    ax.set_xticks(RUN_INDICES, ("Run 1", "Run 2", "Run 3"))
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=4)
    ax.set_yticks(range(4), DATASETS)
    ax.tick_params(axis="y", length=0)
    ax.set_title("Top-ranked item differs across runs", loc="left", pad=8)
    ax.spines["left"].set_visible(False)
    ax.spines["bottom"].set_visible(False)
    handles = [
        Line2D(
            [],
            [],
            marker=SIGN_MARKERS[sign],
            linestyle="none",
            markerfacecolor=SIGN_COLORS[sign],
            markeredgecolor="white",
            label=label,
        )
        for sign, label in ((1, "positive sign"), (-1, "negative sign"))
    ]
    ax.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.06),
        ncol=2,
        borderaxespad=0,
        fontsize=styles.FS_ANNOT - 0.2,
    )


def _panel_b(ax: Axes, data: Figure4Data) -> None:
    x_by_n = {n_value: index for index, n_value in enumerate(N_VALUES)}
    dataset_offsets = np.linspace(-0.24, 0.24, len(DATASETS))
    pair_offsets = {(1, 2): -0.025, (1, 3): 0.0, (2, 3): 0.025}
    for dataset, dataset_offset in zip(DATASETS, dataset_offsets, strict=True):
        subset = data.jaccard[data.jaccard["dataset"] == dataset]
        for row in subset.itertuples(index=False):
            pair = (int(row.run_a_index), int(row.run_b_index))
            ax.scatter(
                x_by_n[int(row.N)] + dataset_offset + pair_offsets[pair],
                float(row.jaccard),
                s=19,
                marker=PAIR_MARKERS[pair],
                facecolor=styles.serotype_color(dataset),
                edgecolor="white",
                linewidth=0.35,
                alpha=0.9,
                zorder=3,
            )
    for y in (0.0, 0.25, 0.5, 0.75):
        ax.axhline(y, color="#E2E2E2", lw=0.45, zorder=0)
    ax.set_xlim(-0.55, 2.55)
    ax.set_ylim(-0.025, 0.75)
    ax.set_xticks(range(3), [str(value) for value in N_VALUES])
    ax.set_yticks([0.0, 0.25, 0.5, 0.75])
    ax.set_xlabel("Top-N set size")
    ax.set_ylabel("Pairwise top-N Jaccard overlap")
    ax.set_title("Within-dataset top-N set overlap", loc="left")

    dataset_handles = [
        Line2D(
            [],
            [],
            marker="o",
            linestyle="none",
            markerfacecolor=styles.serotype_color(dataset),
            markeredgecolor="white",
            label=dataset,
        )
        for dataset in DATASETS
    ]
    pair_handles = [
        Line2D(
            [],
            [],
            marker=marker,
            linestyle="none",
            markerfacecolor="#777777",
            markeredgecolor="white",
            label=f"{a}-{b}",
        )
        for (a, b), marker in PAIR_MARKERS.items()
    ]
    first = ax.legend(
        handles=dataset_handles,
        loc="upper left",
        ncol=2,
        title="Dataset",
        title_fontsize=styles.FS_ANNOT,
        fontsize=styles.FS_ANNOT - 0.3,
    )
    ax.add_artist(first)
    ax.legend(
        handles=pair_handles,
        loc="upper right",
        ncol=1,
        title="Run pair",
        title_fontsize=styles.FS_ANNOT,
        fontsize=styles.FS_ANNOT - 0.3,
    )


def _panel_c(ax: Axes, data: Figure4Data, dataset: str) -> None:
    ranks, signs, labels = heatmap_matrix(data, dataset)
    cmap = LinearSegmentedColormap.from_list(
        "figure4_rank", ["#315E6B", "#A8C1C5", "#EDF1F1"], N=10
    ).with_extremes(bad="white")
    norm = BoundaryNorm(np.arange(0.5, 11.5, 1.0), cmap.N)
    ax.imshow(
        np.ma.masked_invalid(ranks),
        cmap=cmap,
        norm=norm,
        aspect="auto",
        interpolation="nearest",
    )
    for row in range(ranks.shape[0] + 1):
        ax.axhline(row - 0.5, color="#E5E5E5", lw=0.25, zorder=2)
    for column in range(4):
        ax.axvline(column - 0.5, color="#D6D6D6", lw=0.3, zorder=2)
    for row, column in np.argwhere(np.isfinite(ranks)):
        rank = int(ranks[row, column])
        sign = int(signs[row, column])
        sign_text = "+" if sign > 0 else "-" if sign < 0 else "0"
        ax.text(
            column,
            row,
            f"{rank}{sign_text}",
            ha="center",
            va="center",
            fontsize=4.6,
            color="white" if rank <= 4 else "#222222",
            zorder=3,
        )
    ax.set_xticks(range(3), ("Run 1", "Run 2", "Run 3"))
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, labelsize=styles.FS_ANNOT - 0.6, pad=2)
    ax.set_yticks(range(len(labels)), labels)
    ax.tick_params(axis="y", length=0, labelsize=4.6, pad=1)
    ax.set_title(
        f"{dataset}  (union n={len(labels)})",
        loc="left",
        fontsize=styles.FS_LABEL,
        pad=7,
    )
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#BDBDBD")
        spine.set_linewidth(0.35)


def make_figure(data: Figure4Data) -> Figure:
    """Construct the locked three-panel Figure 4 from validated Phase 3E data."""

    styles.apply_style()
    fig = plt.figure(figsize=(DOUBLE_COL, 6.8), constrained_layout=False)
    outer = fig.add_gridspec(
        2,
        2,
        height_ratios=(0.78, 1.82),
        left=0.08,
        right=0.985,
        bottom=0.07,
        top=0.91,
        hspace=0.36,
        wspace=0.30,
    )
    ax_a = fig.add_subplot(outer[0, 0])
    ax_b = fig.add_subplot(outer[0, 1])
    _panel_a(ax_a, data)
    _panel_b(ax_b, data)
    panel_letter(ax_a, "A", dx=-0.09, dy=1.02)
    panel_letter(ax_b, "B", dx=-0.09, dy=1.02)

    lower = outer[1, :].subgridspec(1, 4, wspace=0.55)
    heat_axes: list[Axes] = []
    for index, dataset in enumerate(DATASETS):
        ax = fig.add_subplot(lower[0, index])
        _panel_c(ax, data, dataset)
        heat_axes.append(ax)
    panel_letter(heat_axes[0], "C", dx=-0.42, dy=1.13)
    heat_axes[0].text(
        0.0,
        1.13,
        "Union of top-10 items across runs",
        transform=heat_axes[0].transAxes,
        ha="left",
        va="bottom",
        fontsize=styles.FS_TITLE,
        fontweight="normal",
    )
    heat_axes[-1].text(
        1.0,
        1.13,
        "Cell = rank and sign; blank = outside top 10",
        transform=heat_axes[-1].transAxes,
        ha="right",
        va="bottom",
        fontsize=styles.FS_ANNOT - 0.3,
        color="#666666",
    )

    fig.suptitle(
        "Single-run top-item rankings vary across replicate runs",
        x=0.08,
        y=0.985,
        ha="left",
        fontsize=styles.FS_TITLE,
        fontweight="bold",
    )
    return fig


def numeric_audit(data: Figure4Data) -> dict[str, object]:
    """Return deterministic verification metadata without new inference."""

    top1_rows = []
    for row in data.top1.itertuples(index=False):
        top1_rows.append(
            {
                "dataset": str(row.dataset),
                "run": str(row.run),
                "run_index": int(row.run_index),
                "canonical_coordinate": str(row.canonical_coordinate),
                "effect": float(row.effect),
                "sign": int(row.sign),
                "absolute_effect": float(row.absolute_effect),
            }
        )
    overlap_rows = []
    for row in data.jaccard.itertuples(index=False):
        overlap_rows.append(
            {
                "dataset": str(row.dataset),
                "N": int(row.N),
                "run_a_index": int(row.run_a_index),
                "run_b_index": int(row.run_b_index),
                "intersection_count": int(row.intersection_count),
                "union_count": int(row.union_count),
                "jaccard": float(row.jaccard),
            }
        )
    global_ranges: dict[str, dict[str, float]] = {}
    for n_value in N_VALUES:
        values = data.jaccard.loc[
            data.jaccard["N"] == n_value, "jaccard"
        ].to_numpy(dtype=float)
        global_ranges[str(n_value)] = {
            "min": float(np.min(values)),
            "median": float(np.median(values)),
            "max": float(np.max(values)),
        }
    union_sizes = {
        dataset: int(
            data.heatmap.loc[
                data.heatmap["dataset"] == dataset, "canonical_coordinate"
            ].nunique()
        )
        for dataset in DATASETS
    }
    repeated = (
        data.top1[
            data.top1.duplicated(
                ["dataset", "canonical_coordinate"], keep=False
            )
        ][["dataset", "run_index", "canonical_coordinate"]]
        .to_dict(orient="records")
    )
    return {
        "schema_version": "phase3h.figure4.audit.1",
        "scientific_inputs": list(INPUT_FILES),
        "input_sha256": dict(data.input_hashes),
        "top1_rows": top1_rows,
        "unique_top1_nominations": top1_unique_counts(data),
        "repeated_top1_records": repeated,
        "topn_jaccard_rows": overlap_rows,
        "jaccard_ranges_by_dataset_and_N": jaccard_ranges(data),
        "jaccard_ranges_global_by_N": global_ranges,
        "top10_union_size": union_sizes,
        "tie_count": data.tie_count,
        "rows_in_exact_ties": data.tied_rows,
        "deterministic_tie_rule": data.tie_rule,
        "ranking_basis": "absolute_effect",
        "sign_retained_separately": True,
        "rho_threshold_used": False,
    }


def _caption_text(data: Figure4Data) -> str:
    return (
        "# Figure 4 caption draft\n\n"
        "**Figure 4 | Single-run top-item rankings vary across replicate runs.** "
        "(A) For each frozen dataset and run, the item with the largest absolute "
        "per-run effect is shown; marker orientation retains the effect sign "
        "separately. (B) Within each dataset, Jaccard overlap is shown for all "
        "three run pairs at top-N set sizes of 5, 10, and 20. (C) For the union "
        "of items ranked in the top 10 by at least one run, cells report rank and "
        "sign; blank cells indicate that an item was outside the top 10 in that "
        "run. The figure describes ranking variation in the frozen runs. It does "
        "not identify a true item or provide an inferential estimate of ranking "
        "reliability.\n"
    )


def _write_metadata(paths: Paths, data: Figure4Data) -> list[Path]:
    audit_path = paths.output / "figure4_numeric_audit.json"
    provenance_path = paths.output / "figure4_source_provenance.json"
    caption_path = paths.output / "figure4_caption_draft.md"
    audit_path.write_text(
        json.dumps(numeric_audit(data), indent=2, sort_keys=True) + "\n"
    )
    provenance_path.write_text(
        json.dumps(
            {
                "schema_version": "phase3h.figure4.provenance.1",
                "phase3e_manifest_sha256": PHASE3E_MANIFEST_SHA256,
                "inputs": dict(data.input_hashes),
                "scientific_inputs": list(INPUT_FILES),
                "plotting_module": Path(__file__).name,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n"
    )
    caption_path.write_text(_caption_text(data))
    return [audit_path, provenance_path, caption_path]


def build(paths: Paths) -> list[Path]:
    """Validate Phase 3E inputs, render Figure 4, and write audit metadata."""

    data = load_phase3e_data()
    figure = make_figure(data)
    written = save_figure(figure, paths, FIGURE_STEM)
    _write_metadata(paths, data)
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for path in build(resolve_paths()):
        print(path)
