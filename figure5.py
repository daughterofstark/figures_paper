"""Figure 5: presence and provisional item-level recurrence vary across datasets.

Scientific inputs are restricted to the deterministic Phase 3E presentation
tables in ``derived_data/phase3e``. The renderer verifies their hashes and locked
semantics before displaying:

A. each coordinate-dataset state at the provisional rho=0.5 reference;
B. joint coordinate counts by number of datasets present and crossing;
C. sensitivity of strict all-four operational recurrence to four provisional
   references.

Absence, present-below-reference, and crossing are explicit source states. No
raw frozen table is read and no statistical estimator is recomputed.
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
from matplotlib.colors import BoundaryNorm, LinearSegmentedColormap, ListedColormap
from matplotlib.figure import Figure
from matplotlib.patches import Patch, Rectangle

import styles
from figure_config import DOUBLE_COL, SEROTYPE_ORDER, Paths
from utils import panel_letter, save_figure

FIGURE_STEM = "figure5_cross_serotype"
PHASE3E_MANIFEST_SHA256 = (
    "2043b2c418c5524964492ee32ea9cfad82393b2b8dd60db7b98f1690a6f89873"
)
INPUT_FILES = (
    "figure5_status_matrix.csv",
    "figure5_coordinate_summary.csv",
    "figure5_row_order.csv",
)
EXPECTED_COORDINATES = 248
EXPECTED_PARTIAL_PRESENCE = 49
EXPECTED_ALL_WHERE_PRESENT = 30
EXPECTED_ALL_FOUR = {0.4: 53, 0.5: 18, 0.5330838393: 13, 0.6: 7}
EXPECTED_CROSSING_VECTORS = {
    0.4: (5, 26, 73, 91, 53),
    0.5: (18, 49, 100, 63, 18),
    0.5330838393: (26, 58, 98, 53, 13),
    0.6: (42, 84, 82, 33, 7),
}
DATASETS = tuple(SEROTYPE_ORDER)
STATES = ("ABSENT", "PRESENT_BELOW", "CROSSING")
STATE_CODES = {state: index for index, state in enumerate(STATES)}
STATE_COLORS = {
    "ABSENT": "#FFFFFF",
    "PRESENT_BELOW": "#C9CDD1",
    "CROSSING": "#2A8C82",
}
PROVISIONAL_REFERENCE = 0.5
FORBIDDEN_FIGURE_TERMS = (
    "conserved",
    "conservation",
    "reproducible_all",
    "divergent",
    "shared core",
)

REQUIRED_COLUMNS = {
    "figure5_status_matrix.csv": {
        "cutoff",
        "canonical_coordinate",
        "dataset",
        "chain",
        "configured_group_label",
        "rho",
        "state",
    },
    "figure5_coordinate_summary.csv": {
        "cutoff",
        "canonical_coordinate",
        "chain",
        "configured_group_label",
        "n_total",
        "n_present",
        "n_absent",
        "n_crossing",
        "n_present_below_cutoff",
        "present_in_all_four",
        "crossing_all_four",
    },
    "figure5_row_order.csv": {
        "row_index",
        "canonical_coordinate",
        "chain",
        "configured_group_label",
        "n_present",
        "n_crossing",
    },
}


@dataclass(frozen=True)
class Figure5Data:
    """Validated Phase 3E tables used by Figure 5."""

    status: pd.DataFrame
    summary: pd.DataFrame
    row_order: pd.DataFrame
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


def _at_cutoff(frame: pd.DataFrame, cutoff: float) -> pd.DataFrame:
    return frame[np.isclose(frame["cutoff"], cutoff, rtol=0.0, atol=1e-12)].copy()


def load_phase3e_data(data_dir: Path | None = None) -> Figure5Data:
    """Load and validate only the three locked Phase 3E Figure 5 tables."""

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
            raise FileNotFoundError(f"missing Phase 3E Figure 5 input: {path}")
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

    data = Figure5Data(
        status=frames["figure5_status_matrix.csv"],
        summary=frames["figure5_coordinate_summary.csv"],
        row_order=frames["figure5_row_order.csv"],
        input_hashes=hashes,
    )
    validate_figure5_data(data)
    return data


def validate_figure5_data(data: Figure5Data) -> None:
    """Fail loudly if the registered presence/recurrence semantics drift."""

    order = data.row_order
    if len(order) != EXPECTED_COORDINATES:
        raise ValueError(
            f"expected {EXPECTED_COORDINATES} row-order entries, found {len(order)}"
        )
    if order["canonical_coordinate"].duplicated().any():
        raise ValueError("duplicate coordinate in supplied row order")
    if list(order["row_index"]) != list(range(EXPECTED_COORDINATES)):
        raise ValueError("row_index is not the supplied contiguous 0..247 order")

    status = data.status
    if set(status["state"].astype(str)) != set(STATES):
        raise ValueError(f"unexpected explicit states: {sorted(status['state'].unique())}")
    cutoffs = sorted(status["cutoff"].astype(float).unique())
    if len(cutoffs) != 4 or not np.allclose(
        cutoffs, tuple(EXPECTED_ALL_FOUR), rtol=0.0, atol=1e-12
    ):
        raise ValueError(f"unexpected provisional references: {cutoffs}")
    if status.duplicated(["cutoff", "canonical_coordinate", "dataset"]).any():
        raise ValueError("duplicate cutoff-coordinate-dataset status row")

    for cutoff in EXPECTED_ALL_FOUR:
        subset = _at_cutoff(status, cutoff)
        if len(subset) != EXPECTED_COORDINATES * len(DATASETS):
            raise ValueError(
                f"cutoff {cutoff} does not contain 248 x 4 explicit state rows"
            )
        if set(subset["dataset"].astype(str)) != set(DATASETS):
            raise ValueError(f"cutoff {cutoff} does not contain all four datasets")
        per_coordinate = subset.groupby("canonical_coordinate", sort=False).size()
        if len(per_coordinate) != EXPECTED_COORDINATES or not per_coordinate.eq(4).all():
            raise ValueError(f"cutoff {cutoff} lacks one state per coordinate-dataset")

        absent = subset["state"].eq("ABSENT")
        below = subset["state"].eq("PRESENT_BELOW")
        crossing = subset["state"].eq("CROSSING")
        if not subset.loc[absent, "rho"].isna().all():
            raise ValueError("an explicit ABSENT row has a finite rho")
        if subset.loc[~absent, "rho"].isna().any():
            raise ValueError("a present explicit state has missing rho")
        if not (
            subset.loc[below, "rho"].to_numpy(dtype=float) < float(cutoff)
        ).all():
            raise ValueError("PRESENT_BELOW row does not satisfy its stated reference")
        if not (
            subset.loc[crossing, "rho"].to_numpy(dtype=float) >= float(cutoff)
        ).all():
            raise ValueError("CROSSING row does not satisfy its stated reference")

    summary = data.summary
    if summary.duplicated(["cutoff", "canonical_coordinate"]).any():
        raise ValueError("duplicate cutoff-coordinate summary row")
    for cutoff, expected_vector in EXPECTED_CROSSING_VECTORS.items():
        subset = _at_cutoff(summary, cutoff)
        if len(subset) != EXPECTED_COORDINATES:
            raise ValueError(f"cutoff {cutoff} does not contain 248 summary rows")
        if not subset["n_total"].eq(4).all():
            raise ValueError("n_total is not 4 for every coordinate")
        if not (
            (subset["n_absent"] == 4 - subset["n_present"])
            & (
                subset["n_present_below_cutoff"]
                == subset["n_present"] - subset["n_crossing"]
            )
            & (subset["n_crossing"] <= subset["n_present"])
        ).all():
            raise ValueError("presence/crossing component arithmetic is inconsistent")
        if not (
            subset["crossing_all_four"].astype(bool)
            == ((subset["n_present"] == 4) & (subset["n_crossing"] == 4))
        ).all():
            raise ValueError("crossing_all_four does not implement the strict rule")
        vector = tuple(
            subset["n_crossing"]
            .value_counts()
            .reindex(range(5), fill_value=0)
            .astype(int)
        )
        if vector != expected_vector:
            raise ValueError(
                f"cutoff {cutoff} crossing vector {vector} != {expected_vector}"
            )
        all_four = int(subset["crossing_all_four"].astype(bool).sum())
        if all_four != EXPECTED_ALL_FOUR[cutoff]:
            raise ValueError(
                f"cutoff {cutoff} all-four count {all_four} "
                f"!= {EXPECTED_ALL_FOUR[cutoff]}"
            )

    at_half = _at_cutoff(summary, PROVISIONAL_REFERENCE)
    partial = int((at_half["n_present"] < 4).sum())
    if partial != EXPECTED_PARTIAL_PRESENCE:
        raise ValueError(
            f"partial-presence count {partial} != {EXPECTED_PARTIAL_PRESENCE}"
        )
    all_where_present = int(
        (at_half["n_crossing"] == at_half["n_present"]).sum()
    )
    if all_where_present != EXPECTED_ALL_WHERE_PRESENT:
        raise ValueError(
            f"all-where-present count {all_where_present} "
            f"!= {EXPECTED_ALL_WHERE_PRESENT}"
        )

    supplied = order.set_index("canonical_coordinate")
    half = at_half.set_index("canonical_coordinate")
    if set(supplied.index) != set(half.index):
        raise ValueError("row-order and 0.5 summary coordinate sets differ")
    aligned = half.reindex(supplied.index)
    if not np.array_equal(
        aligned["n_present"].to_numpy(), supplied["n_present"].to_numpy()
    ):
        raise ValueError("row-order n_present does not match 0.5 summary")
    if not np.array_equal(
        aligned["n_crossing"].to_numpy(), supplied["n_crossing"].to_numpy()
    ):
        raise ValueError("row-order n_crossing does not match 0.5 summary")

    half_status = _at_cutoff(status, PROVISIONAL_REFERENCE)
    state_counts = (
        half_status.groupby("canonical_coordinate", sort=False)["state"]
        .value_counts()
        .unstack(fill_value=0)
        .reindex(columns=STATES, fill_value=0)
        .reindex(supplied.index)
    )
    if not np.array_equal(
        4 - state_counts["ABSENT"].to_numpy(), supplied["n_present"].to_numpy()
    ):
        raise ValueError("explicit ABSENT states disagree with n_present")
    if not np.array_equal(
        state_counts["CROSSING"].to_numpy(), supplied["n_crossing"].to_numpy()
    ):
        raise ValueError("explicit CROSSING states disagree with n_crossing")


def status_matrix_at_reference(
    data: Figure5Data, cutoff: float = PROVISIONAL_REFERENCE
) -> np.ndarray:
    """Return the explicit-state matrix in the supplied row and dataset order."""

    subset = _at_cutoff(data.status, cutoff).set_index(
        ["canonical_coordinate", "dataset"]
    )
    keys = [
        (coordinate, dataset)
        for coordinate in data.row_order["canonical_coordinate"].astype(str)
        for dataset in DATASETS
    ]
    states = subset.reindex(keys)["state"]
    if states.isna().any():
        raise ValueError("supplied row order produced a missing matrix state")
    codes = states.map(STATE_CODES)
    return codes.to_numpy(dtype=int).reshape(EXPECTED_COORDINATES, len(DATASETS))


def summary_at_reference(
    data: Figure5Data, cutoff: float = PROVISIONAL_REFERENCE
) -> pd.DataFrame:
    """Return coordinate summaries in the supplied row order."""

    subset = _at_cutoff(data.summary, cutoff).set_index("canonical_coordinate")
    return subset.reindex(data.row_order["canonical_coordinate"]).reset_index()


def presence_count_grid(
    data: Figure5Data, cutoff: float = PROVISIONAL_REFERENCE
) -> np.ndarray:
    """Return the 4 x 5 valid-cell count grid for display."""

    summary = _at_cutoff(data.summary, cutoff)
    grid = np.full((4, 5), np.nan)
    counts = summary.groupby(["n_present", "n_crossing"], sort=True).size()
    for n_present in range(1, 5):
        for n_crossing in range(n_present + 1):
            grid[n_present - 1, n_crossing] = int(
                counts.get((n_present, n_crossing), 0)
            )
    return grid


def all_four_sensitivity(data: Figure5Data) -> pd.DataFrame:
    """Return strict all-four counts at the four registered references."""

    rows = []
    for cutoff in EXPECTED_ALL_FOUR:
        subset = _at_cutoff(data.summary, cutoff)
        rows.append(
            {
                "cutoff": cutoff,
                "count": int(subset["crossing_all_four"].astype(bool).sum()),
            }
        )
    return pd.DataFrame(rows)


def _panel_a(
    fig: Figure,
    spec,
    data: Figure5Data,
) -> tuple[Axes, Axes, Axes]:
    inner = spec.subgridspec(1, 3, width_ratios=(1.0, 0.34, 0.34), wspace=0.16)
    ax = fig.add_subplot(inner[0, 0])
    ax_present = fig.add_subplot(inner[0, 1], sharey=ax)
    ax_crossing = fig.add_subplot(inner[0, 2], sharey=ax)

    matrix = status_matrix_at_reference(data)
    cmap = ListedColormap([STATE_COLORS[state] for state in STATES])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)
    ax.imshow(matrix, cmap=cmap, norm=norm, aspect="auto", interpolation="nearest")
    for row, column in np.argwhere(matrix == STATE_CODES["ABSENT"]):
        ax.add_patch(
            Rectangle(
                (column - 0.5, row - 0.5),
                1,
                1,
                facecolor="none",
                edgecolor="#777777",
                linewidth=0.15,
                hatch="///",
            )
        )
    ax.set_xticks(range(len(DATASETS)), DATASETS, rotation=35, ha="right")
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=2)
    ax.set_yticks([])
    ax.set_ylabel(f"{EXPECTED_COORDINATES} coordinates")
    ax.set_title(r"Status at provisional $\rho=0.5$", loc="left", pad=8)
    for spine in ax.spines.values():
        spine.set_visible(True)
        spine.set_color("#888888")
        spine.set_linewidth(0.4)

    summary = summary_at_reference(data)
    y = np.arange(EXPECTED_COORDINATES)
    ax_present.barh(
        y,
        summary["n_present"],
        height=1.0,
        color="#8D9398",
        edgecolor="none",
    )
    ax_crossing.barh(
        y,
        summary["n_crossing"],
        height=1.0,
        color=STATE_COLORS["CROSSING"],
        edgecolor="none",
    )
    for side_ax, title in (
        (ax_present, "n present"),
        (ax_crossing, "n crossing"),
    ):
        side_ax.set_xlim(0, 4)
        side_ax.set_ylim(EXPECTED_COORDINATES - 0.5, -0.5)
        side_ax.set_xticks([0, 4])
        side_ax.xaxis.tick_top()
        side_ax.tick_params(axis="x", length=2, pad=2)
        side_ax.tick_params(axis="y", left=False, labelleft=False)
        side_ax.set_title(title, fontsize=styles.FS_ANNOT, pad=8)
        side_ax.spines["left"].set_visible(False)
        side_ax.spines["bottom"].set_visible(False)
        side_ax.spines["right"].set_visible(False)

    legend_handles = [
        Patch(
            facecolor=STATE_COLORS["ABSENT"],
            edgecolor="#777777",
            hatch="///",
            linewidth=0.4,
            label="Absent",
        ),
        Patch(
            facecolor=STATE_COLORS["PRESENT_BELOW"],
            edgecolor="none",
            label="Present below reference",
        ),
        Patch(
            facecolor=STATE_COLORS["CROSSING"],
            edgecolor="none",
            label="Crossing",
        ),
    ]
    ax.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(0.0, -0.015),
        fontsize=styles.FS_ANNOT - 0.3,
        ncol=1,
        handlelength=1.2,
        borderaxespad=0,
    )
    panel_letter(ax, "A", dx=-0.10, dy=1.04)
    return ax, ax_present, ax_crossing


def _panel_b(ax: Axes, data: Figure5Data) -> None:
    grid = presence_count_grid(data)
    cmap = LinearSegmentedColormap.from_list(
        "figure5_counts", ["#F1F3F2", "#B9D2CE", "#5C9E98"]
    )
    masked = np.ma.masked_invalid(grid)
    ax.imshow(
        masked,
        origin="lower",
        extent=(-0.5, 4.5, 0.5, 4.5),
        aspect="auto",
        interpolation="nearest",
        cmap=cmap,
        vmin=0,
        vmax=float(np.nanmax(grid)),
    )
    for n_present in range(1, 5):
        for n_crossing in range(n_present + 1):
            value = int(grid[n_present - 1, n_crossing])
            color = "white" if value >= 55 else "#222222"
            ax.text(
                n_crossing,
                n_present,
                str(value),
                ha="center",
                va="center",
                fontsize=styles.FS_LABEL,
                color=color,
                fontweight="bold" if (n_present, n_crossing) == (4, 4) else "normal",
            )

    ax.add_patch(
        Rectangle(
            (3.5, 3.5),
            1,
            1,
            fill=False,
            edgecolor="#111111",
            linewidth=1.1,
            clip_on=False,
        )
    )
    ax.text(
        3.15,
        1.25,
        "Outlined cell: all-four operational recurrence\n"
        "Everywhere present = 30; all four = 18",
        ha="center",
        va="center",
        fontsize=styles.FS_ANNOT - 0.2,
        color="#555555",
    )
    ax.set_xlim(-0.5, 4.5)
    ax.set_ylim(0.5, 4.5)
    ax.set_xticks(range(5))
    ax.set_yticks(range(1, 5))
    ax.set_xlabel(r"$n_{\mathrm{crossing}}$")
    ax.set_ylabel(r"$n_{\mathrm{present}}$")
    ax.set_title(r"Presence-aware counts at provisional $\rho=0.5$", loc="left")
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def _panel_c(ax: Axes, data: Figure5Data) -> None:
    sensitivity = all_four_sensitivity(data)
    x = sensitivity["cutoff"].to_numpy(dtype=float)
    y = sensitivity["count"].to_numpy(dtype=float)
    ax.plot(
        x,
        y,
        color="#555555",
        lw=0.9,
        marker="o",
        markersize=4.2,
        markerfacecolor=STATE_COLORS["CROSSING"],
        markeredgecolor=STATE_COLORS["CROSSING"],
    )
    offsets = ((0, 5), (0, 5), (0, -12), (0, 5))
    for xv, yv, offset in zip(x, y, offsets, strict=True):
        ax.annotate(
            str(int(yv)),
            (xv, yv),
            xytext=offset,
            textcoords="offset points",
            ha="center",
            va="bottom" if offset[1] >= 0 else "top",
            fontsize=styles.FS_LABEL,
            fontweight="bold",
            color="#222222",
        )
    ax.set_xlim(0.385, 0.615)
    ax.set_ylim(0, 60)
    ax.set_xticks(x, ["0.400", "0.500", "0.533", "0.600"], rotation=30, ha="right")
    ax.set_yticks([0, 20, 40, 60])
    ax.set_xlabel(r"Provisional $\rho$ reference")
    ax.set_ylabel("Coordinates")
    ax.set_title("All-four operational recurrence sensitivity", loc="left")


def make_figure(data: Figure5Data) -> Figure:
    """Construct the locked three-panel Figure 5 from validated data."""

    styles.apply_style()
    fig = plt.figure(figsize=(DOUBLE_COL, 5.5), constrained_layout=False)
    outer = fig.add_gridspec(
        1,
        2,
        width_ratios=(1.08, 1.0),
        left=0.075,
        right=0.985,
        bottom=0.13,
        top=0.82,
        wspace=0.34,
    )
    _panel_a(fig, outer[0, 0], data)

    right = outer[0, 1].subgridspec(2, 1, height_ratios=(1.18, 0.88), hspace=0.48)
    ax_b = fig.add_subplot(right[0, 0])
    ax_c = fig.add_subplot(right[1, 0])
    _panel_b(ax_b, data)
    _panel_c(ax_c, data)
    panel_letter(ax_b, "B", dx=-0.10, dy=1.03)
    panel_letter(ax_c, "C", dx=-0.10, dy=1.03)

    fig.suptitle(
        "Presence and provisional item-level recurrence vary across datasets",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=styles.FS_TITLE,
        fontweight="bold",
    )
    return fig


def numeric_audit(data: Figure5Data) -> dict[str, object]:
    """Return deterministic Figure 5 verification metadata."""

    half_summary = _at_cutoff(data.summary, PROVISIONAL_REFERENCE)
    half_status = _at_cutoff(data.status, PROVISIONAL_REFERENCE)

    counts_by_present = {
        str(index): int(value)
        for index, value in half_summary["n_present"]
        .value_counts()
        .reindex(range(1, 5), fill_value=0)
        .sort_index()
        .items()
    }
    counts_by_crossing = {
        str(index): int(value)
        for index, value in half_summary["n_crossing"]
        .value_counts()
        .reindex(range(5), fill_value=0)
        .sort_index()
        .items()
    }
    status_totals = {
        state: int((half_status["state"] == state).sum()) for state in STATES
    }
    all_where_present = int(
        (half_summary["n_crossing"] == half_summary["n_present"]).sum()
    )
    crossing_all_four = int(half_summary["crossing_all_four"].astype(bool).sum())

    all_four_by_cutoff: list[dict[str, object]] = []
    vectors: list[dict[str, object]] = []
    for cutoff in EXPECTED_ALL_FOUR:
        subset = _at_cutoff(data.summary, cutoff)
        vector = [
            int(value)
            for value in subset["n_crossing"]
            .value_counts()
            .reindex(range(5), fill_value=0)
            .sort_index()
        ]
        all_four_by_cutoff.append(
            {
                "cutoff": cutoff,
                "count": int(subset["crossing_all_four"].astype(bool).sum()),
            }
        )
        vectors.append({"cutoff": cutoff, "counts_for_0_to_4_crossings": vector})

    return {
        "schema_version": "phase3g.figure5.audit.1",
        "scientific_inputs": list(INPUT_FILES),
        "input_sha256": dict(sorted(data.input_hashes.items())),
        "total_coordinates": int(half_summary["canonical_coordinate"].nunique()),
        "partial_presence_coordinates": int((half_summary["n_present"] < 4).sum()),
        "counts_by_n_present": counts_by_present,
        "counts_by_n_crossing_at_0.5": counts_by_crossing,
        "status_state_totals_at_0.5": status_totals,
        "crossing_all_four_at_0.5": crossing_all_four,
        "crossing_everywhere_present_at_0.5": all_where_present,
        "difference_all_where_present_minus_all_four": (
            all_where_present - crossing_all_four
        ),
        "all_four_by_cutoff": all_four_by_cutoff,
        "crossing_count_vectors": vectors,
    }


def _write_json(path: Path, payload: object) -> None:
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def write_metadata(data: Figure5Data, output_dir: Path) -> list[Path]:
    """Write deterministic Figure 5 audit, provenance, and caption artifacts."""

    audit = numeric_audit(data)
    audit_path = output_dir / "figure5_numeric_audit.json"
    _write_json(audit_path, audit)

    provenance_path = output_dir / "figure5_source_provenance.json"
    _write_json(
        provenance_path,
        {
            "schema_version": "phase3g.figure5.provenance.1",
            "phase3e_manifest_sha256": PHASE3E_MANIFEST_SHA256,
            "scientific_inputs": [
                {"path": filename, "sha256": data.input_hashes[filename]}
                for filename in INPUT_FILES
            ],
        },
    )

    caption = (
        "# Figure 5 caption draft\n\n"
        "**Figure 5 | Presence and provisional item-level recurrence vary across "
        "datasets.** Here, `n_total=4`; `n_present` is the number of datasets in "
        "which a coordinate exists; `n_absent=4-n_present`; `n_crossing` is the "
        "number of present datasets in which item rho meets the stated provisional "
        "reference; and `n_present_below_cutoff=n_present-n_crossing`. "
        "`crossing_all_four` requires `n_present=4` and `n_crossing=4`. "
        "(A) Explicit coordinate-dataset states at provisional rho=0.5 distinguish "
        "absence from presence below the reference and from crossing; side strips "
        "show the corresponding presence and crossing counts. (B) Joint counts by "
        "`n_present` and `n_crossing` show that crossing in every dataset where a "
        f"coordinate is present ({audit['crossing_everywhere_present_at_0.5']}) is "
        "not equivalent to being present and crossing in all four "
        f"({audit['crossing_all_four_at_0.5']}). (C) The strict all-four operational "
        "set changes across four predeclared provisional references. Absence is not "
        "treated as a below-reference observation. V9B validated deterministic "
        "counting and field semantics only; it did not validate homology or "
        "biological conservation. The cross-dataset coordinate mapping lacks an "
        "archived independent alignment artifact.\n"
    )
    caption_path = output_dir / "figure5_caption_draft.md"
    caption_path.write_text(caption)
    return [audit_path, provenance_path, caption_path]


def build(
    paths: Paths,
    *,
    data_dir: Path | None = None,
    write_supporting_metadata: bool = True,
) -> list[Path]:
    """Render Figure 5 from verified Phase 3E inputs."""

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
