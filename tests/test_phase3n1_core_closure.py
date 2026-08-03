"""Regression checks for the Phase 3N.1 stress-test core closure."""
from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
FROZEN = ROOT / "strideouts"
MANUSCRIPT = ROOT / "STRIDE_MANUSCRIPT_DESCRIPTIVE_METHODS_WORKING.md"


def _bh(values: list[float]) -> list[float]:
    indexed = [(i, p) for i, p in enumerate(values) if np.isfinite(p)]
    result = [float("nan")] * len(values)
    running = float("inf")
    m = len(indexed)
    indexed.sort(key=lambda pair: pair[1])
    for rank in range(m, 0, -1):
        index, p_value = indexed[rank - 1]
        running = min(running, p_value * m / rank)
        result[index] = min(1.0, max(0.0, running))
    return result


def _wald_p(beta: float, beta_se: float) -> float:
    if not np.isfinite(beta) or not np.isfinite(beta_se) or beta_se <= 0:
        return float("nan")
    z_value = abs(beta / beta_se)
    cdf = 0.5 * (1.0 + math.erf(z_value / math.sqrt(2.0)))
    return 2.0 * (1.0 - cdf)


def _selected_rows() -> pd.DataFrame:
    stride = pd.read_parquet(FROZEN / "outputs" / "stride_table.parquet")
    return stride[stride["is_gated_scale"].astype(bool)].copy()


def test_alias_collapse_is_structural_exact_support_deduplication() -> None:
    stride = pd.read_parquet(FROZEN / "outputs" / "stride_table.parquet")
    named_total = 0
    unique_total = 0
    static_duplicate_rho_checks = 0
    for _dataset, frame in stride.groupby("serotype", sort=True):
        support_by_object: dict[tuple[str, str], frozenset[str]] = {}
        rho_by_object: dict[tuple[str, str], float] = {}
        for key, rows in frame.groupby(["scale_level", "region_id"], sort=False):
            support_by_object[key] = frozenset(rows["canon_label"].astype(str))
            rho_values = rows["rho"].drop_duplicates().to_numpy(float)
            assert len(rho_values) == 1
            rho_by_object[key] = float(rho_values[0])

        named_total += len(support_by_object)
        unique_supports = set(support_by_object.values())
        unique_total += len(unique_supports)
        for support in unique_supports:
            objects = [
                key for key, members in support_by_object.items()
                if members == support
            ]
            reference = rho_by_object[objects[0]]
            for key in objects[1:]:
                assert abs(rho_by_object[key] - reference) <= 1e-12
                static_duplicate_rho_checks += 1

        for _item, rows in frame.groupby("canon_label", sort=False):
            ordered = rows.sort_values("scale_index", kind="stable")
            named = [
                (
                    support_by_object[(row.scale_level, row.region_id)],
                    float(row.rho),
                )
                for row in ordered.itertuples(index=False)
            ]
            collapsed = []
            seen: set[frozenset[str]] = set()
            for support, rho in named:
                if support not in seen:
                    seen.add(support)
                    collapsed.append((support, rho))
            named_crossing = next(
                (row for row in named if np.isfinite(row[1]) and row[1] >= 0.5),
                None,
            )
            collapsed_crossing = next(
                (
                    row for row in collapsed
                    if np.isfinite(row[1]) and row[1] >= 0.5
                ),
                None,
            )
            assert named_crossing == collapsed_crossing

    assert (named_total, unique_total, static_duplicate_rho_checks) == (
        1036,
        964,
        72,
    )


def test_interval_exclusion_rows_map_to_174_distinct_fits() -> None:
    selected = _selected_rows()[
        ["serotype", "canon_label", "mech_region_id"]
    ].drop_duplicates(["serotype", "canon_label"])
    screen = pd.read_parquet(
        FROZEN / "outputs_s4" / "significance_screen.parquet"
    )
    joined = screen.merge(
        selected,
        on=["serotype", "canon_label"],
        validate="one_to_one",
    )
    mask = joined["ci_excludes_zero"].astype(bool)
    fit_id = (
        joined.loc[mask, "serotype"].astype(str)
        + "::"
        + joined.loc[mask, "mech_region_id"].astype(str)
    )
    assert int(mask.sum()) == 181
    assert fit_id.nunique() == 174


def test_eight_decimal_rounding_changes_no_frozen_bh_flag() -> None:
    selected = _selected_rows()[
        [
            "serotype",
            "canon_label",
            "mech_direction",
            "mech_beta_signed",
            "mech_beta_se",
        ]
    ].drop_duplicates(["serotype", "canon_label"])
    screen = pd.read_parquet(
        FROZEN / "outputs_s4" / "significance_screen.parquet"
    )
    rows = screen.merge(
        selected,
        on=["serotype", "canon_label"],
        validate="one_to_one",
    )
    rows["p_unrounded"] = [
        _wald_p(beta, beta_se) if direction != "mixed" else float("nan")
        for beta, beta_se, direction in zip(
            rows["mech_beta_signed"],
            rows["mech_beta_se"],
            rows["mech_direction"],
            strict=True,
        )
    ]
    rows["p_rounded"] = rows["p_unrounded"].map(
        lambda value: round(float(value), 8) if np.isfinite(value) else value
    )
    np.testing.assert_allclose(
        rows["p_rounded"],
        rows["p_value"],
        rtol=0,
        atol=0,
        equal_nan=True,
    )

    rows["bh_implemented"] = np.nan
    rows["bh_unrounded"] = np.nan
    for _dataset, indices in rows.groupby("serotype", sort=False).groups.items():
        idx = list(indices)
        implemented = _bh(rows.loc[idx, "p_rounded"].tolist())
        unrounded = _bh(rows.loc[idx, "p_unrounded"].tolist())
        rows.loc[idx, "bh_implemented"] = [
            round(value, 8) if np.isfinite(value) else value
            for value in implemented
        ]
        rows.loc[idx, "bh_unrounded"] = unrounded
    np.testing.assert_allclose(
        rows["bh_implemented"],
        rows["p_value_bh"],
        rtol=0,
        atol=0,
        equal_nan=True,
    )
    implemented_flags = rows["is_signed"] & rows["bh_implemented"].le(0.05)
    unrounded_flags = rows["is_signed"] & rows["bh_unrounded"].le(0.05)
    assert implemented_flags.sum() == unrounded_flags.sum() == 152
    assert np.array_equal(implemented_flags, unrounded_flags)


def test_core_prose_distinguishes_structural_and_empirical_alias_evidence() -> None:
    text = MANUSCRIPT.read_text()
    assert "72,000 checks across the 4,000 V9A dataset realizations" in text
    assert "not Monte Carlo invariance rates" in text
    assert "conditions on registered alias equivalence classes" in text
    assert "consistency guards, not independent statistical validation" in text
    assert "0 of 924,000 decisions" not in text
