"""Focused tests for the Phase 3E presentation-only preparation layer."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import prepare_phase3e_data as p3e


@pytest.fixture(scope="session")
def authoritative_sources() -> p3e.SourcePaths:
    frozen = Path(__file__).resolve().parents[2] / "strideouts"
    return p3e.default_sources(frozen)


@pytest.fixture(scope="session")
def phase3e_build(
    tmp_path_factory: pytest.TempPathFactory,
    authoritative_sources: p3e.SourcePaths,
) -> tuple[Path, dict[str, Path]]:
    output = tmp_path_factory.mktemp("phase3e") / "derived"
    written = p3e.prepare_all(authoritative_sources, output)
    return output, written


def _read(output: Path, name: str) -> pd.DataFrame:
    return pd.read_csv(output / name, float_precision="round_trip")


def test_unique_supports_reproduce_profiles_and_aliases(
    phase3e_build: tuple[Path, dict[str, Path]],
) -> None:
    output, _ = phase3e_build
    unique = _read(output, "figure2_unique_supports.csv")
    assert unique.loc[unique["is_singleton"], ["dataset", "canonical_item_key"]].shape[0] == 924
    assert unique.loc[
        unique["is_first_unique_aggregate"], ["dataset", "canonical_item_key"]
    ].shape[0] == 924
    assert len(unique) == 3493
    patterns = unique["original_named_levels"].value_counts().to_dict()
    assert patterns == {
        "residue": 924,
        "protein|complex": 924,
        "secondary_structure|motif|domain": 721,
        "chain": 721,
        "secondary_structure|motif|domain|chain": 203,
    }
    assert not unique["collapsed_presentation_level"].isin(
        ["secondary_structure", "motif", "domain", "protein", "complex"]
    ).any()


def test_collapse_preserves_every_named_rho_and_membership(
    phase3e_build: tuple[Path, dict[str, Path]],
    authoritative_sources: p3e.SourcePaths,
) -> None:
    output, _ = phase3e_build
    unique = _read(output, "figure2_unique_supports.csv")
    stride = pd.read_parquet(authoritative_sources.stride_table)
    rho_lookup = {
        (str(row.serotype), str(row.canon_label), str(row.scale_level)): float(row.rho)
        for row in stride.itertuples(index=False)
    }
    for row in unique.itertuples(index=False):
        levels = str(row.original_named_levels).split("|")
        source_values = [
            rho_lookup[(str(row.dataset), str(row.canonical_item_key), level)]
            for level in levels
        ]
        assert all(value == source_values[0] for value in source_values)
        assert float(row.rho) == source_values[0]
        members = str(row.support_membership).split(";")
        assert p3e.membership_hash(members) == row.support_membership_hash
        assert len(members) == int(row.support_size)


def test_first_unique_aggregate_is_deterministic(
    phase3e_build: tuple[Path, dict[str, Path]],
) -> None:
    output, _ = phase3e_build
    paired = _read(output, "figure2_item_first_aggregate.csv")
    assert len(paired) == 924
    assert paired.duplicated(["dataset", "canonical_item_key"]).sum() == 0
    assert (paired["first_aggregate_support_size"] > 1).all()
    np.testing.assert_array_equal(
        paired["delta_rho"].to_numpy(),
        (
            paired["first_aggregate_rho"].to_numpy()
            - paired["item_rho"].to_numpy()
        ),
    )


def test_v9a_tables_trace_exactly_to_completed_artifacts(
    authoritative_sources: p3e.SourcePaths,
) -> None:
    support, association, aggregate, crossing, alias = p3e.prepare_v9a_tables(
        authoritative_sources
    )
    results = json.loads(authoritative_sources.v9a_results.read_text())
    source_association = {
        (p3e.DATASET_ID_MAP[row["dataset_id"]], row["regime"]): row
        for row in results["support_size_association"]
    }
    for row in association.itertuples(index=False):
        source = source_association[(row.dataset, row.regime)]
        assert row.mean_spearman == source["mean_coefficient"]
        assert row.mc_t95_lower == source["mc_t95_lower"]
        assert row.mc_t95_upper == source["mc_t95_upper"]
    assert len(support) == 48
    assert len(crossing) == 40
    assert dict(zip(alias["metric"], alias["count"], strict=True)) == {
        "total_decisions": 924000,
        "labels_changed": 672379,
        "selected_membership_changed": 0,
        "selected_rho_changed_gt_1e-12": 0,
    }
    zero_mean = aggregate[aggregate["regime"] == "zero_mean"]
    assert len(zero_mean) == 4
    assert (zero_mean["fraction_increased_mean"] > 0.9999).all()
    assert set(crossing["crossing_level"]) == {
        "item",
        "configured_group",
        "distinct_chain",
        "root",
        "unresolved",
    }


def test_figure5_registered_counts_and_states(
    phase3e_build: tuple[Path, dict[str, Path]],
) -> None:
    output, _ = phase3e_build
    status = _read(output, "figure5_status_matrix.csv")
    summary = _read(output, "figure5_coordinate_summary.csv")
    assert status["canonical_coordinate"].nunique() == 248
    assert len(status) == 248 * 4 * 4
    assert set(status["state"]) == {"ABSENT", "PRESENT_BELOW", "CROSSING"}
    assert int((summary[np.isclose(summary["cutoff"], 0.5)]["n_present"] < 4).sum()) == 49
    for cutoff, expected in p3e.EXPECTED_FIGURE5_VECTORS.items():
        subset = summary[np.isclose(summary["cutoff"], cutoff)]
        vector = tuple(
            subset["n_crossing"]
            .value_counts()
            .reindex(range(5), fill_value=0)
            .astype(int)
        )
        assert vector == expected
        assert int(subset["crossing_all_four"].sum()) == expected[4]
        assert (
            subset["crossing_all_four"]
            == ((subset["n_present"] == 4) & (subset["n_crossing"] == 4))
        ).all()
    assert (
        status["state"].eq("ABSENT")
        == status["rho"].isna()
    ).all()
    present = status[status["state"] != "ABSENT"]
    assert (
        present["state"].eq("CROSSING")
        == (present["rho"] >= present["cutoff"])
    ).all()


def test_figure5_row_order_is_explicit_and_stable(
    phase3e_build: tuple[Path, dict[str, Path]],
) -> None:
    output, _ = phase3e_build
    order = _read(output, "figure5_row_order.csv")
    assert list(order["row_index"]) == list(range(248))
    expected = order.copy()
    expected["_number"] = expected["canonical_coordinate"].map(
        lambda value: p3e.canonical_sort_key(value)[1]
    )
    expected = expected.sort_values(
        [
            "n_present",
            "n_crossing",
            "chain",
            "_number",
            "canonical_coordinate",
        ],
        ascending=[False, False, True, True, True],
        kind="stable",
    )
    assert list(expected["canonical_coordinate"]) == list(order["canonical_coordinate"])


def test_figure4_ranking_and_jaccard_are_descriptive_only(
    phase3e_build: tuple[Path, dict[str, Path]],
) -> None:
    output, _ = phase3e_build
    top1 = _read(output, "figure4_top1.csv")
    jaccard = _read(output, "figure4_topn_jaccard.csv")
    heatmap = _read(output, "figure4_top10_heatmap.csv")
    assert len(top1) == 12
    assert set(jaccard["N"]) == {5, 10, 20}
    assert (jaccard.groupby(["dataset", "N"]).size() == 3).all()
    row = jaccard.iloc[0]
    assert row["jaccard"] == row["intersection_count"] / row["union_count"]
    assert (heatmap.groupby(["dataset", "canonical_coordinate"]).size() == 3).all()
    assert (
        heatmap["in_top10"].astype(bool)
        == (heatmap["rank"] <= 10)
    ).all()
    forbidden = {"rho", "cutoff", "supportable_scale", "recovered"}
    for frame in (top1, jaccard, heatmap):
        assert forbidden.isdisjoint(frame.columns)


def test_real_ranking_has_no_exact_absolute_effect_ties(
    phase3e_build: tuple[Path, dict[str, Path]],
) -> None:
    output, _ = phase3e_build
    manifest = json.loads((output / "manifest.json").read_text())
    assert manifest["ranking"]["exact_tie_groups"] == 0
    assert manifest["ranking"]["rows_in_exact_ties"] == 0
    assert manifest["ranking"]["tie_break"] == "canonical_coordinate_ascending"


def test_tie_break_fixture_uses_canonical_key() -> None:
    fixture = pd.DataFrame(
        {
            "serotype": ["D1"] * 3,
            "replicate": ["r1"] * 3,
            "replicate_index": [1] * 3,
            "canon_label": ["X:2", "X:1", "X:3"],
            "r": [-0.5, 0.5, 0.2],
            "chain": ["X"] * 3,
            "domain": ["G"] * 3,
        }
    )
    ranked, diagnostics = p3e._rank_effects(fixture, Path("fixture"))
    assert list(ranked["canonical_coordinate"])[:2] == ["X:1", "X:2"]
    assert list(ranked["rank"]) == [1, 2, 3]
    assert diagnostics == {"exact_tie_groups": 1, "rows_in_exact_ties": 2}
    assert list(ranked["sign"])[:2] == [1, -1]


def test_configured_group_support_sizes_are_exact(
    phase3e_build: tuple[Path, dict[str, Path]],
    authoritative_sources: p3e.SourcePaths,
) -> None:
    output, _ = phase3e_build
    derived = _read(output, "configured_group_supports.csv")
    source = pd.read_parquet(authoritative_sources.replicate_table)
    expected = (
        source[["serotype", "canon_label", "chain", "domain", "region_id"]]
        .drop_duplicates(["serotype", "canon_label"])
        .groupby(["serotype", "chain", "domain", "region_id"])
        .size()
        .rename("support_size")
        .reset_index()
        .sort_values(["serotype", "chain", "domain"], kind="stable")
        .reset_index(drop=True)
    )
    assert len(derived) == 32
    assert list(derived["support_size"]) == list(expected["support_size"])
    assert list(derived["configured_group_label"]) == list(expected["domain"])


def test_same_inputs_produce_byte_identical_outputs(
    tmp_path: Path,
    authoritative_sources: p3e.SourcePaths,
) -> None:
    first = tmp_path / "first"
    second = tmp_path / "second"
    one = p3e.prepare_all(authoritative_sources, first)
    two = p3e.prepare_all(authoritative_sources, second)
    assert set(one) == set(two)
    for filename in sorted(one):
        assert one[filename].read_bytes() == two[filename].read_bytes(), filename
