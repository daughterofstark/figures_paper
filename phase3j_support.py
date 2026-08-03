"""Shared Phase 3J presentation metadata and deterministic audit helpers."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import pandas as pd

_HERE = Path(__file__).resolve().parent
_PHASE3E_SUPPORTS = _HERE / "derived_data" / "phase3e" / "configured_group_supports.csv"
_PHASE3E_MANIFEST = _HERE / "derived_data" / "phase3e" / "manifest.json"


def sha256(path: Path) -> str:
    """Return the SHA-256 digest of *path*."""
    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def write_json(path: Path, payload: object) -> None:
    """Write deterministic, human-readable JSON."""
    path.write_text(json.dumps(payload, indent=2, sort_keys=True, allow_nan=False) + "\n")


def load_configured_group_metadata(
    observed: pd.DataFrame,
) -> tuple[pd.DataFrame, dict[str, str]]:
    """Load locked Phase 3E metadata when keys match the authoritative dataset.

    The synthetic regression fixture has a different key set. For that fixture
    only, support sizes are deterministically counted from its item-level S7
    table by the calling script and passed through ``observed``.
    """
    required = {"serotype", "chain", "domain", "support_size"}
    missing = required - set(observed.columns)
    if missing:
        raise ValueError(f"support metadata input missing columns: {sorted(missing)}")

    observed_keys = set(
        observed[["serotype", "chain", "domain"]]
        .astype(str)
        .itertuples(index=False, name=None)
    )
    if _PHASE3E_SUPPORTS.is_file() and _PHASE3E_MANIFEST.is_file():
        manifest = json.loads(_PHASE3E_MANIFEST.read_text())
        expected = manifest["outputs"]["configured_group_supports.csv"]["sha256"]
        actual = sha256(_PHASE3E_SUPPORTS)
        if actual != expected:
            raise ValueError(
                "Phase 3E configured-group metadata hash mismatch: "
                f"expected {expected}, observed {actual}"
            )
        locked = pd.read_csv(_PHASE3E_SUPPORTS)
        locked_keys = set(
            locked[["dataset", "chain", "configured_group_label"]]
            .astype(str)
            .itertuples(index=False, name=None)
        )
        if observed_keys == locked_keys:
            out = locked.rename(
                columns={
                    "dataset": "serotype",
                    "configured_group_label": "domain",
                }
            ).copy()
            return out, {
                "metadata_source": str(_PHASE3E_SUPPORTS),
                "metadata_sha256": actual,
            }

    # Synthetic/generic test path: the supplied support_size is already an
    # exact deterministic count from the fixture's item membership.
    out = observed.drop_duplicates(["serotype", "chain", "domain"]).copy()
    out["is_configured_catalytic_group"] = out.get(
        "is_configured_catalytic_group", False
    )
    return out, {
        "metadata_source": "input item membership (generic fixture fallback)",
        "metadata_sha256": "",
    }


def group_support_counts(item_rows: pd.DataFrame) -> pd.DataFrame:
    """Count exact item rows per dataset/configured-group key."""
    required = {"serotype", "chain", "domain", "canon_label"}
    missing = required - set(item_rows.columns)
    if missing:
        raise ValueError(f"item table missing support columns: {sorted(missing)}")
    if item_rows.duplicated(["serotype", "canon_label"]).any():
        raise ValueError("item table has duplicate dataset-item rows")
    return (
        item_rows.groupby(["serotype", "chain", "domain"], sort=True)["canon_label"]
        .nunique()
        .rename("support_size")
        .reset_index()
    )
