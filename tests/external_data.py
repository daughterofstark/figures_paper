"""Shared discovery and skip policy for archived-data integration tests."""
from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path

import pytest

from figure_config import ExternalSources, resolve_external_sources


def sources() -> ExternalSources:
    return resolve_external_sources()


def require_files(paths: Iterable[Path], purpose: str) -> None:
    missing = [str(path) for path in paths if not path.is_file()]
    if missing:
        pytest.skip(
            f"integration data unavailable for {purpose}; set the documented "
            f"STRIDE_* source variables; missing: {missing}"
        )


def require_frozen_root(root: Path, purpose: str) -> None:
    require_files(
        [
            root / "outputs" / "stride_table.parquet",
            root / "outputs" / "replicate_table.parquet",
        ],
        purpose,
    )
