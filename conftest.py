"""Pytest bootstrap: put the project root (and tests/) on sys.path.

The figure modules use flat imports (``import styles``), so the project directory
must be importable when the tests run.
"""
from __future__ import annotations

import os
import sys
from pathlib import Path

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.join(_HERE, "tests")):
    if _p not in sys.path:
        sys.path.insert(0, _p)


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Skip archived-data integration tests uniformly when inputs are absent."""
    from figure_config import resolve_external_sources

    sources = resolve_external_sources()
    required: tuple[Path, ...] = (
        sources.outputs_root / "outputs" / "stride_table.parquet",
        sources.outputs_root / "outputs" / "replicate_table.parquet",
        sources.manuscript_path,
        sources.v9a_artifacts / "v9a_frozen_design_results.json",
        sources.v9b_artifacts / "v9b_partial_presence_results.json",
    )
    missing = [str(path) for path in required if not path.is_file()]
    if not missing:
        return
    reason = (
        "archived STRIDE integration inputs unavailable; set the documented "
        "STRIDE_* source variables; missing: " + ", ".join(missing)
    )
    marker = pytest.mark.skip(reason=reason)
    for item in items:
        if "integration" in item.keywords:
            item.add_marker(marker)
