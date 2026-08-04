from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import supplementary_figure2 as s2
from figure_config import resolve_paths
from external_data import require_frozen_root, sources

FROZEN = sources().outputs_root
pytestmark = pytest.mark.integration


@pytest.fixture(scope="module")
def signed_data() -> s2.SignedFigureData:
    require_frozen_root(FROZEN, "Supplementary Figure S2")
    paths = resolve_paths(root=FROZEN)
    return s2.load_signed_figure_data(paths)


def test_authoritative_counts_reproduce_exactly(
    signed_data: s2.SignedFigureData,
) -> None:
    assert signed_data.authoritative
    assert s2.bookkeeping_counts(signed_data) == {
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


def test_distinct_fit_identity_uses_explicit_selected_support(
    signed_data: s2.SignedFigureData,
) -> None:
    rows = signed_data.selected_rows
    expected = rows["serotype"].astype(str) + "::" + rows["mech_region_id"].astype(str)
    assert np.array_equal(rows["fit_id"], expected)
    assert signed_data.distinct_selected_fits["fit_id"].is_unique
    assert len(signed_data.distinct_selected_fits) == 538
    assert (
        rows.groupby("fit_id")["beta_signed"].nunique(dropna=False) <= 1
    ).all()
    assert (
        rows.groupby("fit_id")["beta_se"].nunique(dropna=False) <= 1
    ).all()


def test_repeated_representations_are_explicit_and_exact(
    signed_data: s2.SignedFigureData,
) -> None:
    rows = signed_data.selected_rows
    finite = rows[rows["finite_signed_test"]]
    repeated_fits = finite.groupby("fit_id").size()
    assert int((repeated_fits > 1).sum()) == 7
    assert int((repeated_fits - 1).clip(lower=0).sum()) == 34
    assert rows["repeated_representation"].equals(
        rows["item_row_multiplicity"] > 1
    )
    matrix, repeated = s2.coordinate_status_matrix(signed_data)
    assert matrix.shape == repeated.shape == (4, 248)
    assert int(repeated.sum()) == int(
        (rows["repeated_representation"] & rows["finite_signed_test"]).sum()
    )


def test_raw_interval_and_bh_counts_are_not_conflated(
    signed_data: s2.SignedFigureData,
) -> None:
    counts = s2.bookkeeping_counts(signed_data)
    assert counts["raw_interval_exclusion_rows"] == 181
    assert counts["distinct_raw_interval_exclusion_fits"] == 174
    assert counts["raw_wald_flagged_rows"] == 181
    assert counts["bh_flagged_rows"] == 152
    assert counts["raw_interval_exclusion_rows"] != counts["bh_flagged_rows"]


def test_dataset_level_bh_labels_and_families_reproduce(
    signed_data: s2.SignedFigureData,
) -> None:
    rows = signed_data.selected_rows
    finite = rows[rows["finite_signed_test"]]
    assert finite.groupby("serotype").size().to_dict() == {
        "DENV1": 70,
        "DENV2": 79,
        "DENV3": 94,
        "DENV4": 59,
    }
    assert np.array_equal(
        finite["bh_flagged"].to_numpy(),
        (finite["p_value_bh"] <= 0.05).to_numpy(),
    )
    assert (
        signed_data.distinct_selected_fits[
            signed_data.distinct_selected_fits["bh_flagged"]
        ]
        .groupby("serotype")
        .size()
        .to_dict()
        == {"DENV1": 40, "DENV2": 12, "DENV3": 48, "DENV4": 45}
    )


def test_sign_beta_and_intervals_are_unchanged(
    signed_data: s2.SignedFigureData,
) -> None:
    rows = signed_data.selected_rows
    finite = rows[rows["finite_signed_test"]]
    assert np.array_equal(
        np.sign(finite["beta_signed"]).astype(int),
        np.where(finite["direction"] == "increase", 1, -1),
    )
    frozen = pd.read_parquet(
        FROZEN / "outputs_s4" / "significance_screen.parquet"
    ).sort_values(["serotype", "chain", "canon_label"], kind="stable")
    for column in (
        "beta_signed",
        "beta_se",
        "beta_ci_lower",
        "beta_ci_upper",
        "p_value",
        "p_value_bh",
        "significant_fdr",
    ):
        left = rows[column].reset_index(drop=True)
        right = frozen[column].reset_index(drop=True)
        assert left.equals(right), column


def test_panel_a_uses_one_mark_per_distinct_finite_fit(
    signed_data: s2.SignedFigureData,
) -> None:
    assert len(signed_data.distinct_finite_fits) == 268
    assert {
        dataset: len(s2.distinct_fits_by_dataset(signed_data, dataset))
        for dataset in s2.DATASETS
    } == {"DENV1": 66, "DENV2": 57, "DENV3": 89, "DENV4": 56}


def test_panel_b_retains_all_rows_and_explicit_absence(
    signed_data: s2.SignedFigureData,
) -> None:
    matrix, _ = s2.coordinate_status_matrix(signed_data)
    assert matrix.shape == (4, 248)
    counts = dict(zip(*np.unique(matrix, return_counts=True), strict=True))
    assert counts[0] == 68
    assert sum(counts.values()) == 4 * 248
    assert sum(value for key, value in counts.items() if key != 0) == 924


def test_unsafe_figure_facing_terminology_is_absent(
    signed_data: s2.SignedFigureData,
) -> None:
    fig = s2.make_figure(signed_data)
    text = " ".join(
        artist.get_text()
        for artist in fig.findobj(
            lambda obj: hasattr(obj, "get_text") and callable(obj.get_text)
        )
    ).lower()
    plt.close(fig)
    for forbidden in s2.FORBIDDEN_FIGURE_TERMS:
        assert forbidden not in text


def test_build_writes_deterministic_figure_metadata_and_caption(
    tmp_path: Path,
) -> None:
    first = resolve_paths(root=FROZEN, output=tmp_path / "first")
    second = resolve_paths(root=FROZEN, output=tmp_path / "second")
    written_first = s2.build(first)
    written_second = s2.build(second)
    assert {path.suffix for path in written_first} == {".svg", ".pdf", ".png"}
    for filename in (
        "Supplementary_Figure_S2.svg",
        "Supplementary_Figure_S2.pdf",
        "Supplementary_Figure_S2.png",
        "Supplementary_Figure_S2_numeric_audit.json",
        "Supplementary_Figure_S2_source_provenance.json",
        "Supplementary_Figure_S2_caption_draft.md",
    ):
        assert (first.output / filename).read_bytes() == (
            second.output / filename
        ).read_bytes()
    audit = json.loads(
        (first.output / "Supplementary_Figure_S2_numeric_audit.json").read_text()
    )
    assert audit["counts"]["selected_item_rows"] == 924
    assert audit["counts"]["distinct_finite_fits"] == 268
    assert audit["post_selection_fdr_control_established"] is False
