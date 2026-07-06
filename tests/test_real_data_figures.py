"""Regression tests for real-data (multi-chain) STRIDE outputs.

These reproduce the exact Figure 6 failure — a domain label that recurs across
chains within a serotype, which made ``domain`` a non-unique index — and assert
the whole project now builds on that real-like data. They also guard Figure 3
(which must not average distinct chains) and a full ``build_all`` end-to-end.
"""
from __future__ import annotations

import warnings
from pathlib import Path

import pandas as pd
import pytest

from realistic_data import write_realistic_outputs

from figure_config import resolve_paths


@pytest.fixture()
def realistic_paths(tmp_path: Path):
    root = write_realistic_outputs(tmp_path / "outs")
    return resolve_paths(root=root, output=tmp_path / "figures")


def test_duplicate_domain_is_present_in_fixture(realistic_paths):
    """Guard: the fixture really does have a domain on two chains within a serotype."""
    df = pd.read_parquet(realistic_paths.s7("F6_variance_composition.parquet"))
    d1 = df[df["serotype"] == "DENV1"]
    dup = d1["domain"].value_counts()
    assert (dup > 1).any(), "fixture must contain a duplicated domain label"
    # ...and it is disambiguated by chain
    assert not d1.duplicated(["chain", "domain"]).any()


def test_figure6_builds_on_duplicate_domain(realistic_paths):
    """The exact previously-crashing case: reindex on duplicate 'domain' labels."""
    import figure6

    written = figure6.build(realistic_paths)
    assert {p.suffix for p in written} == {".svg", ".pdf", ".png"}
    for p in written:
        assert p.is_file() and p.stat().st_size > 0


def test_figure3_builds_and_keeps_both_chains(realistic_paths):
    import figure3

    written = figure3.build(realistic_paths)
    assert len(written) == 3
    # the heatmap must retain one row per (chain, domain), not collapse them
    df = pd.read_parquet(realistic_paths.s7("F3_domain_serotype_rho_heatmap.parquet"))
    regions = df[["chain", "domain"]].drop_duplicates()
    assert (regions["domain"].value_counts() > 1).any()  # a domain spans 2 chains


def test_figure8_builds_on_duplicate_domain(realistic_paths):
    import figure8

    assert len(figure8.build(realistic_paths)) == 3


def test_full_build_all_on_realistic_outputs(realistic_paths):
    """End-to-end: every figure + table builds with no exception and no font warning."""
    import build_all

    with warnings.catch_warnings():
        warnings.simplefilter("error", UserWarning)
        results = build_all.build_all(realistic_paths)

    total = sum(len(v) for v in results.values())
    assert total == 8 * 3 + 5 * 3  # 24 figure files + 15 table files
    # spot-check a representative artifact of each kind exists
    out = realistic_paths.output
    assert (out / "figure6_variance_composition.svg").is_file()
    assert (out / "table5_variance_component_budget.tex").is_file()
