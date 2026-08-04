from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import tables
from figure_config import resolve_paths
from external_data import require_frozen_root, sources

pytestmark = pytest.mark.integration

ROOT = sources().outputs_root
FORBIDDEN = (
    "fdr-controlled",
    "robust effect",
    "significant mechanism",
    "supported scale",
    "licensed scale",
    "conserved position",
    "reproducible_all",
    "biological domain",
    "causal variance",
    "mechanistic effect",
)


def _prepared() -> tables.PreparedTables:
    require_frozen_root(ROOT, "repaired Tables 1-5")
    return tables.prepare_tables(resolve_paths(root=ROOT))


def test_table1_distinguishes_interval_exclusion_and_bh_flags() -> None:
    table1 = _prepared().frames["table1_per_serotype_summary"]
    assert list(table1["dataset"]) == ["DENV1", "DENV2", "DENV3", "DENV4"]
    assert table1["item_rows"].sum() == 924
    assert table1["item_rows_rho_ge_provisional_0.5"].to_list() == [
        122,
        137,
        145,
        106,
    ]
    assert table1["selected_signed_rows_with_finite_Wald_test"].sum() == 302
    assert table1["selected_signed_rows_Wald_interval_excludes_zero"].sum() == 181
    assert _prepared().consistency["distinct_raw_interval_exclusion_fits"] == 174
    assert table1["BH_flagged_selected_item_rows"].sum() == 152
    assert (
        table1["selected_signed_rows_Wald_interval_excludes_zero"].sum()
        != table1["BH_flagged_selected_item_rows"].sum()
    )
    assert "provisional" in tables.TABLE_SPECS[
        "table1_per_serotype_summary"
    ]["note"].lower()


def test_table2_has_exact_groups_supports_and_distinct_fit_bookkeeping() -> None:
    prepared = _prepared()
    table2 = prepared.frames["table2_configured_group_summary"]
    signed = prepared.frames["supplementary_table2_conditional_signed_bookkeeping"]
    assert len(table2) == 32
    assert len(signed) == 32
    assert table2["support_size"].sum() == 924
    assert np.array_equal(table2["support_size"], table2["selected_item_rows"])
    assert table2["support_size"].min() == 3
    assert table2["support_size"].max() == 131
    assert "distinct_selected_fits_represented" in table2
    assert {
        "finite_signed_rows",
        "distinct_finite_fits_represented",
        "BH_flagged_rows",
        "distinct_BH_flagged_fits_represented",
    }.issubset(signed.columns)
    assert prepared.consistency["selected_item_rows"] == 924
    assert prepared.consistency["distinct_selected_fits"] == 538
    assert prepared.consistency["finite_signed_rows"] == 302
    assert prepared.consistency["distinct_finite_fits"] == 268
    assert prepared.consistency["bh_flagged_rows"] == 152
    assert prepared.consistency["distinct_bh_flagged_fits"] == 145


def test_table3_is_minimally_relabelled_with_support_sizes() -> None:
    table3 = _prepared().frames["table3_configured_catalytic_groups"]
    assert len(table3) == 8
    assert set(table3["configured_group"]) == {"Catalytic Triad", "Oxyanion Loop"}
    assert set(table3["support_size"]) == {3, 8}
    assert "magnitude_rho" in table3
    assert "fitted_magnitude_beta" in table3


def test_table4_is_strict_all_four_and_preserves_denominator_semantics() -> None:
    prepared = _prepared()
    table4 = prepared.frames["table4_all_four_operational_recurrence"]
    assert len(table4) == 18
    assert (table4["n_total"] == 4).all()
    assert (table4["n_present"] == 4).all()
    assert (table4["n_absent"] == 0).all()
    assert (table4["n_crossing"] == 4).all()
    assert (table4["n_present_below_cutoff"] == 0).all()
    assert table4["present_in_all_four"].all()
    assert table4["crossing_all_four"].all()
    assert prepared.consistency["union_coordinates"] == 248
    assert prepared.consistency["partial_presence_coordinates"] == 49
    assert prepared.consistency["all_where_present_count"] == 30
    assert prepared.consistency["strict_all_four_count"] == 18


def test_table5_preserves_all_values_supports_and_28_2_2_classes() -> None:
    prepared = _prepared()
    table5 = prepared.frames["table5_configured_group_variance_composition"]
    assert len(table5) == 32
    assert table5["support_size"].between(3, 131).all()
    assert prepared.consistency["variance_regime_counts"] == {
        "balanced": 2,
        "replicate-dominated": 28,
        "sampling-dominated": 2,
    }
    frozen = pd.read_parquet(
        ROOT / "outputs_s7" / "T5_variance_component_budget.parquet"
    ).sort_values(["serotype", "chain", "domain"])
    observed = table5.copy()
    observed[["domain", "chain"]] = observed["configured_group"].str.extract(
        r"^(.*) \(([^()]*)\)$"
    )
    observed = observed.sort_values(["dataset", "chain", "domain"])
    for table_column, frozen_column in (
        ("fitted_between_run_component", "tau2"),
        ("propagated_within_run_component", "sigma2_bar"),
        ("between_run_fraction", "frac_tau2"),
        ("within_run_fraction", "frac_sigma2"),
    ):
        np.testing.assert_array_equal(
            observed[table_column].to_numpy(), frozen[frozen_column].to_numpy()
        )
    note = tables.TABLE_SPECS[
        "table5_configured_group_variance_composition"
    ]["note"].lower()
    assert "k=3" in note
    assert "zero within-support covariance" in note


def test_table_facing_text_uses_safe_terminology() -> None:
    tables._assert_safe_table_text()
    text = " ".join(
        item
        for spec in tables.TABLE_SPECS.values()
        for item in (spec["caption"], spec["note"])
    ).lower()
    for forbidden in FORBIDDEN:
        assert forbidden not in text


def test_repaired_table_exports_are_byte_deterministic(tmp_path: Path) -> None:
    first = resolve_paths(root=ROOT, output=tmp_path / "first")
    second = resolve_paths(root=ROOT, output=tmp_path / "second")
    written_first = tables.build(first)
    written_second = tables.build(second)
    assert len(written_first) == len(written_second) == 24
    for first_path, second_path in zip(written_first, written_second, strict=True):
        assert first_path.name == second_path.name
        assert first_path.read_bytes() == second_path.read_bytes()
    for name in (
        "phase3k_table_manifest.json",
        "phase3k_cross_table_consistency_audit.json",
    ):
        assert (first.output / name).read_bytes() == (second.output / name).read_bytes()
    audit = json.loads(
        (first.output / "phase3k_cross_table_consistency_audit.json").read_text()
    )
    assert audit["raw_interval_exclusion_rows"] == 181
    assert audit["distinct_raw_interval_exclusion_fits"] == 174
    assert audit["bh_flagged_rows"] == 152
    assert audit["strict_all_four_count"] == 18
