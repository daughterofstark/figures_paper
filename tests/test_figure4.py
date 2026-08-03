from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pytest

import figure4
from figure_config import resolve_paths


@pytest.fixture(scope="module")
def figure4_data() -> figure4.Figure4Data:
    return figure4.load_phase3e_data()


def test_only_registered_phase3e_inputs_are_declared() -> None:
    assert figure4.INPUT_FILES == (
        "figure4_top1.csv",
        "figure4_topn_jaccard.csv",
        "figure4_top10_heatmap.csv",
    )


def test_panel_a_reproduces_all_top1_entries(
    figure4_data: figure4.Figure4Data,
) -> None:
    top1 = figure4_data.top1
    assert len(top1) == 12
    assert set(top1["dataset"]) == set(figure4.DATASETS)
    assert not top1.duplicated(["dataset", "run_index"]).any()
    assert top1.groupby("dataset")["run_index"].agg(list).to_dict() == {
        dataset: [1, 2, 3] for dataset in figure4.DATASETS
    }
    np.testing.assert_array_equal(
        top1["sign"].astype(int), np.sign(top1["effect"]).astype(int)
    )
    assert figure4.top1_unique_counts(figure4_data) == {
        "DENV1": 2,
        "DENV2": 3,
        "DENV3": 3,
        "DENV4": 3,
    }


def test_panel_b_reproduces_all_pairwise_jaccard_values(
    figure4_data: figure4.Figure4Data,
) -> None:
    overlap = figure4_data.jaccard
    assert set(overlap["N"]) == {5, 10, 20}
    assert len(overlap) == 36
    assert overlap["jaccard"].between(0, 1).all()
    quotient = overlap["intersection_count"] / overlap["union_count"]
    np.testing.assert_allclose(quotient, overlap["jaccard"], rtol=0, atol=1e-15)
    assert (
        overlap.groupby(["dataset", "N"], sort=False).size().eq(3).all()
    )
    assert set(
        zip(
            overlap["run_a_index"].astype(int),
            overlap["run_b_index"].astype(int),
            strict=True,
        )
    ) == {(1, 2), (1, 3), (2, 3)}


def test_panel_b_has_no_inferential_fields_or_logic(
    figure4_data: figure4.Figure4Data,
) -> None:
    forbidden = {"p_value", "p_value_bh", "significant", "ci_lower", "ci_upper"}
    assert forbidden.isdisjoint(figure4_data.jaccard.columns)
    source = Path(figure4.__file__).read_text().lower()
    assert "ttest" not in source
    assert "bootstrap" not in source


def test_panel_c_preserves_union_order_rank_membership_and_sign(
    figure4_data: figure4.Figure4Data,
) -> None:
    heatmap = figure4_data.heatmap
    assert len(heatmap) == 315
    for dataset in figure4.DATASETS:
        subset = heatmap[heatmap["dataset"] == dataset]
        ranks, signs, labels = figure4.heatmap_matrix(figure4_data, dataset)
        source_labels = (
            subset[["union_order", "canonical_coordinate"]]
            .drop_duplicates()
            .sort_values("union_order", kind="stable")["canonical_coordinate"]
            .astype(str)
            .tolist()
        )
        assert labels == source_labels
        assert ranks.shape == signs.shape == (len(labels), 3)
        for row in subset.itertuples(index=False):
            i = int(row.union_order)
            j = int(row.run_index) - 1
            if bool(row.in_top10):
                assert ranks[i, j] == int(row.rank)
                assert signs[i, j] == int(row.sign)
            else:
                assert np.isnan(ranks[i, j])
                assert signs[i, j] == 0


def test_tie_provenance_is_explicit(
    figure4_data: figure4.Figure4Data,
) -> None:
    assert figure4_data.tie_count == 0
    assert figure4_data.tied_rows == 0
    assert figure4_data.tie_rule == "canonical_coordinate_ascending"


def test_no_rho_or_threshold_logic_contaminates_inputs(
    figure4_data: figure4.Figure4Data,
) -> None:
    for frame in (
        figure4_data.top1,
        figure4_data.jaccard,
        figure4_data.heatmap,
    ):
        assert {
            "rho",
            "cutoff",
            "threshold",
            "p_value",
            "significant",
        }.isdisjoint(frame.columns)
    assert figure4.numeric_audit(figure4_data)["rho_threshold_used"] is False


def test_figure_facing_text_contains_no_retired_terms(
    figure4_data: figure4.Figure4Data,
) -> None:
    fig = figure4.make_figure(figure4_data)
    text = " ".join(
        artist.get_text()
        for artist in fig.findobj(
            lambda obj: hasattr(obj, "get_text") and callable(obj.get_text)
        )
    ).lower()
    plt.close(fig)
    for forbidden in figure4.FORBIDDEN_FIGURE_TERMS:
        assert forbidden not in text


def test_hash_validation_rejects_modified_input(
    tmp_path: Path,
) -> None:
    source = figure4.default_data_dir()
    for name in (*figure4.INPUT_FILES, "manifest.json"):
        (tmp_path / name).write_bytes((source / name).read_bytes())
    target = tmp_path / "figure4_top1.csv"
    target.write_bytes(target.read_bytes() + b"\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        figure4.load_phase3e_data(tmp_path)


def test_build_writes_only_figure4_and_deterministic_metadata(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    paths1 = resolve_paths(root=tmp_path / "unused", output=first)
    paths2 = resolve_paths(root=tmp_path / "unused", output=second)
    written1 = figure4.build(paths1)
    written2 = figure4.build(paths2)
    assert {path.suffix for path in written1} == {".svg", ".pdf", ".png"}
    for name in (
        "figure4_numeric_audit.json",
        "figure4_source_provenance.json",
        "figure4_caption_draft.md",
    ):
        assert (first / name).read_bytes() == (second / name).read_bytes()
    audit = json.loads((first / "figure4_numeric_audit.json").read_text())
    assert audit["scientific_inputs"] == list(figure4.INPUT_FILES)
    assert audit["unique_top1_nominations"] == {
        "DENV1": 2,
        "DENV2": 3,
        "DENV3": 3,
        "DENV4": 3,
    }
    assert not (first / "figure2_spatial_resolution.png").exists()
    assert not (first / "figure5_cross_serotype.png").exists()
