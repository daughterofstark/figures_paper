from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import figure5
from figure_config import resolve_paths


@pytest.fixture(scope="module")
def figure5_data() -> figure5.Figure5Data:
    return figure5.load_phase3e_data()


def test_only_registered_phase3e_inputs_are_declared() -> None:
    assert figure5.INPUT_FILES == (
        "figure5_status_matrix.csv",
        "figure5_coordinate_summary.csv",
        "figure5_row_order.csv",
    )


def test_panel_a_retains_explicit_three_state_matrix_and_supplied_order(
    figure5_data: figure5.Figure5Data,
) -> None:
    order = figure5_data.row_order
    assert len(order) == 248
    assert list(order["row_index"]) == list(range(248))
    matrix = figure5.status_matrix_at_reference(figure5_data)
    assert matrix.shape == (248, 4)
    assert set(np.unique(matrix)) == set(figure5.STATE_CODES.values())
    half = figure5._at_cutoff(figure5_data.status, 0.5)
    assert len(half) == 248 * 4
    assert not half.duplicated(["canonical_coordinate", "dataset"]).any()
    first = str(order.iloc[0]["canonical_coordinate"])
    expected = (
        half[half["canonical_coordinate"] == first]
        .set_index("dataset")
        .reindex(figure5.DATASETS)["state"]
        .map(figure5.STATE_CODES)
        .to_numpy()
    )
    np.testing.assert_array_equal(matrix[0], expected)
    assert set(half["state"]) == {"ABSENT", "PRESENT_BELOW", "CROSSING"}
    assert figure5.STATE_CODES["ABSENT"] != figure5.STATE_CODES["PRESENT_BELOW"]


def test_panel_b_presence_grid_is_valid_and_sums_to_248(
    figure5_data: figure5.Figure5Data,
) -> None:
    grid = figure5.presence_count_grid(figure5_data)
    assert grid.shape == (4, 5)
    assert int(np.nansum(grid)) == 248
    for n_present in range(1, 5):
        assert np.isnan(grid[n_present - 1, n_present + 1 :]).all()
        assert np.isfinite(grid[n_present - 1, : n_present + 1]).all()
    assert int(grid[3, 4]) == 18
    half = figure5._at_cutoff(figure5_data.summary, 0.5)
    assert (half["n_crossing"] <= half["n_present"]).all()
    all_where_present = int((half["n_crossing"] == half["n_present"]).sum())
    all_four = int(half["crossing_all_four"].astype(bool).sum())
    assert all_where_present == 30
    assert all_four == 18
    assert all_where_present != all_four


def test_panel_c_all_four_counts_and_cutoffs_are_exact(
    figure5_data: figure5.Figure5Data,
) -> None:
    sensitivity = figure5.all_four_sensitivity(figure5_data)
    np.testing.assert_allclose(
        sensitivity["cutoff"], [0.4, 0.5, 0.5330838393, 0.6], rtol=0, atol=1e-12
    )
    assert list(sensitivity["count"]) == [53, 18, 13, 7]


def test_all_crossing_count_vectors_are_exact(
    figure5_data: figure5.Figure5Data,
) -> None:
    for cutoff, expected in figure5.EXPECTED_CROSSING_VECTORS.items():
        subset = figure5._at_cutoff(figure5_data.summary, cutoff)
        observed = tuple(
            subset["n_crossing"]
            .value_counts()
            .reindex(range(5), fill_value=0)
            .astype(int)
        )
        assert observed == expected


def test_registered_presence_and_state_counts(
    figure5_data: figure5.Figure5Data,
) -> None:
    audit = figure5.numeric_audit(figure5_data)
    assert audit["total_coordinates"] == 248
    assert audit["partial_presence_coordinates"] == 49
    assert audit["counts_by_n_present"] == {"1": 3, "2": 13, "3": 33, "4": 199}
    assert audit["status_state_totals_at_0.5"] == {
        "ABSENT": 68,
        "PRESENT_BELOW": 414,
        "CROSSING": 510,
    }
    assert audit["difference_all_where_present_minus_all_four"] == 12


def test_explicit_absence_is_not_derived_as_below_cutoff(
    figure5_data: figure5.Figure5Data,
) -> None:
    half = figure5._at_cutoff(figure5_data.status, 0.5)
    absent = half["state"] == "ABSENT"
    below = half["state"] == "PRESENT_BELOW"
    assert absent.any() and below.any()
    assert not (absent & below).any()
    assert half.loc[absent, "rho"].isna().all()
    assert half.loc[below, "rho"].notna().all()
    assert (half.loc[below, "rho"] < 0.5).all()


def test_figure_facing_text_contains_no_retired_terms(
    figure5_data: figure5.Figure5Data,
) -> None:
    fig = figure5.make_figure(figure5_data)
    text = " ".join(
        artist.get_text()
        for artist in fig.findobj(
            lambda obj: hasattr(obj, "get_text") and callable(obj.get_text)
        )
    ).lower()
    plt.close(fig)
    for forbidden in figure5.FORBIDDEN_FIGURE_TERMS:
        assert forbidden not in text


def test_hash_validation_rejects_modified_input(
    tmp_path: Path, figure5_data: figure5.Figure5Data
) -> None:
    source = figure5.default_data_dir()
    for name in (*figure5.INPUT_FILES, "manifest.json"):
        (tmp_path / name).write_bytes((source / name).read_bytes())
    target = tmp_path / "figure5_status_matrix.csv"
    target.write_bytes(target.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        figure5.load_phase3e_data(tmp_path)


def test_build_writes_only_figure5_and_deterministic_metadata(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    paths1 = resolve_paths(root=tmp_path / "unused", output=first)
    paths2 = resolve_paths(root=tmp_path / "unused", output=second)
    written1 = figure5.build(paths1)
    written2 = figure5.build(paths2)
    assert {path.suffix for path in written1} == {".svg", ".pdf", ".png"}
    for name in (
        "figure5_numeric_audit.json",
        "figure5_source_provenance.json",
        "figure5_caption_draft.md",
    ):
        assert (first / name).read_bytes() == (second / name).read_bytes()
    audit = json.loads((first / "figure5_numeric_audit.json").read_text())
    assert audit["total_coordinates"] == 248
    assert audit["scientific_inputs"] == list(figure5.INPUT_FILES)
    assert not (first / "figure2_spatial_resolution.png").exists()
