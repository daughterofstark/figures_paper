from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import figure1
import figure3
import figure6
import supplementary_figure1 as s1
from figure_config import resolve_paths
from external_data import require_frozen_root, sources

pytestmark = pytest.mark.integration

ROOT = sources().outputs_root
FORBIDDEN = (
    "rescue",
    "recovery",
    "emergence",
    "saturation",
    "supported scale",
    "licensed",
    "robust",
    "mechanism",
    "causal",
    "conserved",
    "biological domain",
    "statistically significant",
)


@pytest.fixture(scope="module")
def frozen_paths():
    require_frozen_root(ROOT, "Figures 1, 3, 6, and S1")
    return resolve_paths(root=ROOT)


def test_figure1_preserves_rows_presence_and_across_present_denominator(
    frozen_paths,
) -> None:
    items, presence, order = figure1.prepare_data(frozen_paths)
    assert len(items) == 924
    assert len(order) == 248
    assert items.groupby("serotype").size().to_dict() == {
        "DENV1": 236,
        "DENV2": 247,
        "DENV3": 242,
        "DENV4": 199,
    }
    counts = items.groupby("canon_label")["serotype"].nunique().reindex(order)
    assert np.array_equal(
        counts.to_numpy(),
        presence.set_index("canon_label").reindex(order)["n_serotypes_present"].to_numpy(),
    )
    assert int((counts < 4).sum()) == 49
    assert int((4 - counts).sum()) == 68
    assert not (items["rho_residue"] == 0).all()
    medians = items.groupby("canon_label")["rho_residue"].median().reindex(order)
    np.testing.assert_allclose(
        medians,
        presence.set_index("canon_label").reindex(order)["_rho_across_present"],
        rtol=0,
        atol=0,
    )


def test_s1_profiles_are_exact_and_unsmoothed(frozen_paths) -> None:
    items, order, catalytic = s1.prepare_data(frozen_paths)
    frozen = pd.read_parquet(
        frozen_paths.s7("F1_reproducibility_landscape.parquet")
    )
    assert len(order) == 248
    assert len(catalytic) == 3
    assert items.sort_values(["serotype", "canon_label"]).reset_index(drop=True).equals(
        frozen.sort_values(["serotype", "canon_label"]).reset_index(drop=True)
    )
    source = Path(s1.__file__).read_text()
    assert "_rolling" not in source
    assert ".plot(" not in source
    assert "no interpolation or smoothing" in source
    assert "provisional rho reference = 0.5" in s1.REFERENCE_LABEL


def test_figure3_configured_group_values_and_support_sizes_are_exact(
    frozen_paths,
) -> None:
    data, provenance = figure3.prepare_data(frozen_paths)
    assert len(data) == 32
    assert provenance["metadata_source"].endswith("configured_group_supports.csv")
    assert not data.duplicated(["serotype", "chain", "domain"]).any()
    assert data["support_size"].min() == 3
    assert data["support_size"].max() == 131
    assert int(data["is_configured_catalytic_group"].sum()) == 8
    frozen_heat = pd.read_parquet(
        frozen_paths.s7("F3_domain_serotype_rho_heatmap.parquet")
    ).sort_values(["serotype", "chain", "domain"])
    merged = data.sort_values(["serotype", "chain", "domain"])
    np.testing.assert_array_equal(
        merged["rho_domain"].to_numpy(), frozen_heat["rho_domain"].to_numpy()
    )
    frozen_coherence = pd.read_parquet(
        frozen_paths.s7("F8_coherence_vs_rho.parquet")
    ).sort_values(["serotype", "chain", "domain"])
    np.testing.assert_array_equal(
        merged["coherence_domain"].to_numpy(),
        frozen_coherence["coherence_domain"].to_numpy(),
    )


def test_figure6_values_regimes_and_support_sizes_are_exact(frozen_paths) -> None:
    data, provenance = figure6.prepare_data(frozen_paths)
    assert len(data) == 32
    assert provenance["metadata_source"].endswith("configured_group_supports.csv")
    assert data["variance_regime"].value_counts().to_dict() == {
        "replicate_dominated": 28,
        "balanced": 2,
        "sampling_dominated": 2,
    }
    np.testing.assert_allclose(data["frac_tau2"] + data["frac_sigma2"], 1.0)
    frozen = pd.read_parquet(
        frozen_paths.s7("F6_variance_composition.parquet")
    ).sort_values(["serotype", "chain", "domain"])
    observed = data.sort_values(["serotype", "chain", "domain"])
    for column in ("frac_tau2", "frac_sigma2"):
        np.testing.assert_array_equal(
            observed[column].to_numpy(), frozen[column].to_numpy()
        )
    assert observed["support_size"].between(3, 131).all()


@pytest.mark.parametrize(
    "strings",
    [
        figure1.FIGURE_STRINGS,
        s1.FIGURE_STRINGS,
        figure3.FIGURE_STRINGS,
        figure6.FIGURE_STRINGS,
    ],
)
def test_phase3j_figure_facing_terminology_is_safe(strings) -> None:
    text = " ".join(strings).lower()
    for forbidden in FORBIDDEN:
        assert forbidden not in text


@pytest.mark.parametrize(
    ("module", "stem", "audit", "caption"),
    [
        (
            figure1,
            "figure1_reproducibility_landscape",
            "figure1_numeric_audit.json",
            "figure1_caption_draft.md",
        ),
        (
            s1,
            "Supplementary_Figure_S1",
            "Supplementary_Figure_S1_numeric_audit.json",
            "Supplementary_Figure_S1_caption_draft.md",
        ),
        (
            figure3,
            "figure3_catalytic_core",
            "figure3_numeric_audit.json",
            "figure3_caption_draft.md",
        ),
        (
            figure6,
            "figure6_variance_composition",
            "figure6_numeric_audit.json",
            "figure6_caption_draft.md",
        ),
    ],
)
def test_phase3j_outputs_and_metadata_are_deterministic(
    tmp_path: Path, module, stem: str, audit: str, caption: str
) -> None:
    first = resolve_paths(root=ROOT, output=tmp_path / f"{stem}-first")
    second = resolve_paths(root=ROOT, output=tmp_path / f"{stem}-second")
    module.build(first)
    module.build(second)
    for filename in (
        f"{stem}.svg",
        f"{stem}.pdf",
        f"{stem}.png",
        audit,
        caption,
    ):
        assert (first.output / filename).read_bytes() == (
            second.output / filename
        ).read_bytes()
    payload = json.loads((first.output / audit).read_text())
    assert payload
