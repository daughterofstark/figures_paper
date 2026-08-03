from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pytest

import figure2
from figure_config import resolve_paths


@pytest.fixture(scope="module")
def figure2_data() -> figure2.Figure2Data:
    return figure2.load_phase3e_data()


def test_only_registered_phase3e_inputs_are_declared() -> None:
    assert figure2.INPUT_FILES == (
        "figure2_item_first_aggregate.csv",
        "figure2_unique_supports.csv",
        "figure2_v9a_support_size.csv",
        "figure2_v9a_support_size_association.csv",
        "figure2_v9a_item_first_aggregate.csv",
        "figure2_v9a_first_crossing.csv",
        "figure2_v9a_alias_provenance.csv",
        "figure2_observed_first_crossing_sensitivity.csv",
    )


def test_panel_a_has_exactly_924_unique_observed_profiles(
    figure2_data: figure2.Figure2Data,
) -> None:
    paired = figure2_data.item_aggregate
    assert len(paired) == 924
    assert paired.duplicated(["dataset", "canonical_item_key"]).sum() == 0
    assert set(paired["dataset"]) == set(figure2.EXPECTED_DATASETS)
    assert {"item_rho", "first_aggregate_rho"}.issubset(paired.columns)


def test_panel_b_uses_stored_descriptive_delta_only(
    figure2_data: figure2.Figure2Data,
) -> None:
    paired = figure2_data.item_aggregate
    np.testing.assert_allclose(
        paired["delta_rho"],
        paired["first_aggregate_rho"] - paired["item_rho"],
        rtol=0,
        atol=1e-15,
    )
    assert paired.groupby("dataset").size().to_dict() == {
        "DENV1": 236,
        "DENV2": 247,
        "DENV3": 242,
        "DENV4": 199,
    }
    audit = figure2.numeric_audit(figure2_data)
    assert all("p_value" not in row for row in audit["panel_b"])


def test_panel_c_traces_v9a_values_and_keeps_regimes_separate(
    figure2_data: figure2.Figure2Data,
) -> None:
    support = figure2_data.v9a_support_size
    association = figure2_data.v9a_association
    assert set(support["regime"]) == set(figure2.REGIME_ORDER)
    assert set(association["regime"]) == set(figure2.REGIME_ORDER)
    assert not support.duplicated(
        ["dataset", "regime", "support_size_stratum"]
    ).any()
    assert (
        support["mc_t95_lower"].le(support["mean_rho"])
        & support["mean_rho"].le(support["mc_t95_upper"])
    ).all()
    zero = association[association["regime"] == "zero_mean"].set_index("dataset")
    expected = {
        "DENV1": 0.3347612681735509,
        "DENV2": 0.32786544822679187,
        "DENV3": 0.3305342899653259,
        "DENV4": 0.36039193785456575,
    }
    for dataset, value in expected.items():
        assert np.isclose(
            zero.loc[dataset, "mean_spearman"], value, rtol=0, atol=1e-15
        )
    ordered = support.sort_values(
        ["dataset", "regime", "support_size_stratum_index"], kind="stable"
    )
    for _, rows in ordered.groupby(["dataset", "regime"], sort=False):
        assert rows["support_size_stratum_index"].is_monotonic_increasing


def test_panel_d_has_four_references_and_collapsed_categories_only(
    figure2_data: figure2.Figure2Data,
) -> None:
    crossing = figure2_data.observed_crossing
    assert np.allclose(
        sorted(crossing["cutoff"].unique()), figure2.CUTOFFS, rtol=0, atol=1e-12
    )
    assert set(crossing["crossing_level"]) == set(figure2.CROSSING_CATEGORIES)
    forbidden = {"secondary_structure", "motif", "domain", "protein", "complex"}
    assert forbidden.isdisjoint(set(crossing["crossing_level"]))
    overall = crossing[crossing["dataset"] == "ALL"]
    assert (overall.groupby("cutoff")["count"].sum() == 924).all()
    assert (overall["denominator"] == 924).all()


def test_alias_provenance_is_exact(figure2_data: figure2.Figure2Data) -> None:
    aliases = dict(
        zip(
            figure2_data.alias_provenance["metric"],
            figure2_data.alias_provenance["count"],
            strict=True,
        )
    )
    assert aliases == {
        "total_decisions": 924000,
        "labels_changed": 672379,
        "selected_membership_changed": 0,
        "selected_rho_changed_gt_1e-12": 0,
    }
    caption = (
        Path(figure2.__file__).parent / "output" / "figure2_caption_draft.md"
    ).read_text()
    assert "structurally" in caption
    assert "Monte Carlo invariance rate" in caption
    assert "72,000 duplicate-representation checks" in caption
    assert "0 of 924,000" not in caption


def test_figure_text_contains_no_retired_scientific_terms(
    figure2_data: figure2.Figure2Data,
) -> None:
    fig = figure2.make_figure(figure2_data)
    text = " ".join(
        artist.get_text()
        for artist in fig.findobj(
            lambda obj: hasattr(obj, "get_text") and callable(obj.get_text)
        )
    ).lower()
    plt.close(fig)
    for forbidden in figure2.FORBIDDEN_SCIENTIFIC_TERMS:
        assert forbidden not in text
    assert "secondary structure" not in text
    assert "motif" not in text
    assert "domain" not in text


def test_hash_validation_rejects_modified_input(
    tmp_path: Path, figure2_data: figure2.Figure2Data
) -> None:
    source = figure2.default_data_dir()
    for name in (*figure2.INPUT_FILES, "manifest.json"):
        (tmp_path / name).write_bytes((source / name).read_bytes())
    target = tmp_path / "figure2_item_first_aggregate.csv"
    target.write_bytes(target.read_bytes() + b"\\n")
    with pytest.raises(ValueError, match="hash mismatch"):
        figure2.load_phase3e_data(tmp_path)


def test_build_writes_only_figure2_and_deterministic_metadata(
    tmp_path: Path,
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    paths1 = resolve_paths(root=tmp_path / "unused", output=first)
    paths2 = resolve_paths(root=tmp_path / "unused", output=second)
    written1 = figure2.build(paths1)
    written2 = figure2.build(paths2)
    assert {path.suffix for path in written1} == {".svg", ".pdf", ".png"}
    for name in (
        "figure2_numeric_audit.json",
        "figure2_source_provenance.json",
        "figure2_caption_draft.md",
    ):
        assert (first / name).read_bytes() == (second / name).read_bytes()
    audit = json.loads((first / "figure2_numeric_audit.json").read_text())
    assert audit["panel_a"]["n"] == 924
    assert audit["scientific_inputs"] == list(figure2.INPUT_FILES)
    assert not (first / "figure1_reproducibility_landscape.png").exists()
