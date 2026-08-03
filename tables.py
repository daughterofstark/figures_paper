"""Phase 3K manuscript tables with descriptive, audited semantics.

The builders read frozen S0/S4/S5/S7 tables plus the locked Phase 3E
presentation metadata. They rename, join, and count existing rows
deterministically; no estimator, threshold, p-value, or adjusted value is
recomputed.

CSV outputs preserve source numeric precision. LaTeX, HTML, and Markdown use
the repository's three-significant-digit display convention.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from figure_config import Paths
from phase3j_support import (
    group_support_counts,
    load_configured_group_metadata,
    sha256,
    write_json,
)
from supplementary_figure2 import bookkeeping_counts, load_signed_figure_data
from utils import export_table, load_s7_figure, load_table, parse_canon_label

PROVISIONAL_RHO = 0.5
OPERATIONAL_COHERENCE = 0.6
PHASE3E_DIR = Path(__file__).resolve().parent / "derived_data" / "phase3e"
PHASE3E_MANIFEST = PHASE3E_DIR / "manifest.json"
AUTHORITATIVE_F1_HASH = (
    "08d5852ed86c2d671f52851920643d87a0a3471e515888d92489765a125e731e"
)

UNSAFE_TABLE_TERMS = (
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

TABLE_SPECS = {
    "table1_per_serotype_summary": {
        "caption": "Dataset summary under the provisional item-rho reference.",
        "label": "tab:dataset_summary",
        "note": (
            "Item-rho crossing uses the provisional reference rho=0.5. "
            "Finite signed rows are selected rows with an implemented Wald test. "
            "Interval exclusion is the raw symmetric Wald interval excluding zero; "
            "across all datasets, 181 such item rows represent 174 distinct fitted "
            "effects. It is distinct from the stored within-dataset BH flag. All "
            "displayed counts use the dataset item-row denominator shown."
        ),
    },
    "table2_configured_group_summary": {
        "caption": (
            "Configured-group magnitude stability and operational coherence."
        ),
        "label": "tab:configured_group",
        "note": (
            "Rho is support-specific magnitude stability. Coherence is the "
            "implemented pooled sign-homogeneity statistic. Row counts are selected "
            "item-row representations after provisional rho=0.5 first crossing; "
            "fit counts are unique selected exact-support identities represented "
            "among those rows. A fit spanning multiple configured groups can appear "
            "in more than one group row, so group-level fit counts are not additive. "
            "Finite-test and stored BH-flag bookkeeping is reported separately."
        ),
    },
    "supplementary_table2_conditional_signed_bookkeeping": {
        "caption": (
            "Conditional selected-support signed-test bookkeeping by configured group."
        ),
        "label": "tab:conditional_signed_bookkeeping",
        "note": (
            "Signed fitting follows provisional rho=0.5 first-crossing selection "
            "and operational coherence eligibility at 0.6. Row counts are selected "
            "item-row representations; fit counts are unique selected exact-support "
            "identities represented among those rows. A fit spanning multiple "
            "configured groups can appear in more than one group row, so group-level "
            "fit counts are not additive. Stored BH flags are conditional descriptive "
            "labels, not an established post-selection error-control guarantee."
        ),
    },
    "table3_configured_catalytic_groups": {
        "caption": "Configured catalytic-group magnitude fits across datasets.",
        "label": "tab:configured_catalytic",
        "note": (
            "Rho and beta are support-level magnitude-fit quantities. Exact support "
            "size is shown for each dataset/configured-group row; no signed or "
            "causal interpretation is attached to fitted magnitude beta."
        ),
    },
    "table4_all_four_operational_recurrence": {
        "caption": (
            "Coordinates present and crossing provisional item rho=0.5 in all "
            "four datasets."
        ),
        "label": "tab:all_four_operational",
        "note": (
            "n_present counts datasets in which a coordinate exists; n_absent = "
            "4 - n_present. n_crossing counts present datasets with item rho at or "
            "above the provisional 0.5 reference; n_present_below_cutoff = "
            "n_present - n_crossing. Strict all-four operational recurrence requires "
            "n_present=4 and n_crossing=4. Absence is not treated as below-cutoff "
            "evidence."
        ),
    },
    "table5_configured_group_variance_composition": {
        "caption": (
            "Configured-group fitted variance composition under the K=3 "
            "diagonal-covariance model."
        ),
        "label": "tab:variance_composition",
        "note": (
            "The fitted between-run component and propagated within-run component "
            "are model-based. Propagation assumes fixed support membership and zero "
            "within-support covariance among item estimates. With K=3, heterogeneity "
            "and variance components are weakly identified. Regime labels are "
            "descriptive bookkeeping: between-run fraction >=0.6, "
            "replicate-dominated; <=0.4, sampling-dominated; otherwise balanced."
        ),
    },
}


@dataclass(frozen=True)
class PreparedTables:
    frames: dict[str, pd.DataFrame]
    source_hashes: dict[str, str]
    consistency: dict[str, object]


def _source_hashes(paths: Paths) -> dict[str, str]:
    sources = {
        "outputs_s7/T1_per_serotype_summary.parquet": paths.s7(
            "T1_per_serotype_summary.parquet"
        ),
        "outputs_s7/T2_domain_rho_signed_effect.parquet": paths.s7(
            "T2_domain_rho_signed_effect.parquet"
        ),
        "outputs_s7/T3_catalytic_cross_serotype.parquet": paths.s7(
            "T3_catalytic_cross_serotype.parquet"
        ),
        "outputs_s7/T5_variance_component_budget.parquet": paths.s7(
            "T5_variance_component_budget.parquet"
        ),
        "outputs_s7/F1_reproducibility_landscape.parquet": paths.s7(
            "F1_reproducibility_landscape.parquet"
        ),
        "outputs_s4/significance_screen.parquet": paths.table(
            "s4", "significance_screen.parquet"
        ),
        "outputs/stride_table.parquet": paths.table("s0", "stride_table.parquet"),
        "outputs_s5/position_conservation.parquet": paths.table(
            "s5", "position_conservation.parquet"
        ),
    }
    for path in sources.values():
        if not path.is_file():
            raise FileNotFoundError(f"missing table source: {path}")
    hashes = {name: sha256(path) for name, path in sources.items()}
    if hashes["outputs_s7/F1_reproducibility_landscape.parquet"] == AUTHORITATIVE_F1_HASH:
        manifest = json.loads(PHASE3E_MANIFEST.read_text())
        for name in ("configured_group_supports.csv", "figure5_coordinate_summary.csv"):
            path = PHASE3E_DIR / name
            expected = manifest["outputs"][name]["sha256"]
            actual = sha256(path)
            if actual != expected:
                raise ValueError(
                    f"Phase 3E source hash mismatch for {name}: "
                    f"expected {expected}, observed {actual}"
                )
            hashes[f"derived_data/phase3e/{name}"] = actual
    return hashes


def _group_metadata(paths: Paths) -> tuple[pd.DataFrame, dict[str, str]]:
    item = load_s7_figure(paths, "F1_reproducibility_landscape")
    counts = group_support_counts(item)
    metadata, provenance = load_configured_group_metadata(counts)
    return metadata, provenance


def _table1(paths: Paths, signed_rows: pd.DataFrame) -> pd.DataFrame:
    source = load_s7_figure(paths, "T1_per_serotype_summary").copy()
    required = {
        "serotype",
        "n_loci",
        "n_gated_residue",
        "n_signed",
        "n_mixed",
        "n_signed_significant",
        "frac_mixed",
        "rho_residue_median",
        "rho_residue_q1",
        "rho_residue_q3",
    }
    if missing := required - set(source.columns):
        raise ValueError(f"Table 1 source missing columns: {sorted(missing)}")
    if source["serotype"].duplicated().any():
        raise ValueError("Table 1 source has duplicate dataset rows")
    screen = (
        signed_rows.groupby("serotype", sort=True)
        .agg(
            finite_signed_rows=("finite_signed_test", "sum"),
            interval_excluding_zero_rows=("raw_interval_exclusion", "sum"),
            bh_flagged_rows=("bh_flagged", "sum"),
        )
        .reset_index()
    )
    out = source.merge(screen, on="serotype", how="left", validate="one_to_one")
    if not np.array_equal(out["n_signed"], out["finite_signed_rows"]):
        raise ValueError("Table 1 finite signed counts disagree with S4")
    if not np.array_equal(
        out["n_signed_significant"], out["interval_excluding_zero_rows"]
    ):
        raise ValueError("Table 1 interval-exclusion counts disagree with S4")
    item = load_s7_figure(paths, "F1_reproducibility_landscape").copy()
    item["rho_residue"] = pd.to_numeric(item["rho_residue"], errors="raise")
    crossings = (
        item.assign(crossing=item["rho_residue"] >= PROVISIONAL_RHO)
        .groupby("serotype", sort=True)["crossing"]
        .sum()
    )
    observed_crossings = out.set_index("serotype")["n_gated_residue"]
    if not crossings.reindex(observed_crossings.index).equals(observed_crossings):
        raise ValueError("Table 1 provisional item-rho crossing counts disagree")
    return out[
        [
            "serotype",
            "n_loci",
            "n_gated_residue",
            "finite_signed_rows",
            "n_mixed",
            "frac_mixed",
            "interval_excluding_zero_rows",
            "bh_flagged_rows",
            "rho_residue_median",
            "rho_residue_q1",
            "rho_residue_q3",
        ]
    ].rename(
        columns={
            "serotype": "dataset",
            "n_loci": "item_rows",
            "n_gated_residue": "item_rows_rho_ge_provisional_0.5",
            "n_mixed": "selected_rows_without_finite_signed_test",
            "frac_mixed": "fraction_without_finite_signed_test",
            "finite_signed_rows": "selected_signed_rows_with_finite_Wald_test",
            "interval_excluding_zero_rows": (
                "selected_signed_rows_Wald_interval_excludes_zero"
            ),
            "bh_flagged_rows": "BH_flagged_selected_item_rows",
            "rho_residue_median": "median_item_rho",
            "rho_residue_q1": "item_rho_q1",
            "rho_residue_q3": "item_rho_q3",
        }
    ).sort_values("dataset", kind="stable").reset_index(drop=True)


def _group_count_frame(rows: pd.DataFrame) -> pd.DataFrame:
    keys = ["serotype", "chain", "domain"]
    base = (
        rows.groupby(keys, sort=True)
        .agg(
            selected_item_rows=("canon_label", "size"),
            distinct_selected_fits_represented=("fit_id", "nunique"),
        )
        .reset_index()
    )
    finite = (
        rows[rows["finite_signed_test"]]
        .groupby(keys, sort=True)
        .agg(
            finite_signed_rows=("canon_label", "size"),
            distinct_finite_fits_represented=("fit_id", "nunique"),
        )
        .reset_index()
    )
    flagged = (
        rows[rows["bh_flagged"]]
        .groupby(keys, sort=True)
        .agg(
            bh_flagged_rows=("canon_label", "size"),
            distinct_bh_flagged_fits_represented=("fit_id", "nunique"),
        )
        .reset_index()
    )
    out = base.merge(finite, on=keys, how="left", validate="one_to_one")
    out = out.merge(flagged, on=keys, how="left", validate="one_to_one")
    count_cols = [
        "finite_signed_rows",
        "distinct_finite_fits_represented",
        "bh_flagged_rows",
        "distinct_bh_flagged_fits_represented",
    ]
    out[count_cols] = out[count_cols].fillna(0).astype(int)
    return out


def _table2(
    paths: Paths, signed_rows: pd.DataFrame, metadata: pd.DataFrame
) -> pd.DataFrame:
    source = load_s7_figure(paths, "T2_domain_rho_signed_effect").copy()
    keys = ["serotype", "chain", "domain"]
    required = set(keys) | {"rho_domain", "coherence_domain"}
    if missing := required - set(source.columns):
        raise ValueError(f"Table 2 source missing columns: {sorted(missing)}")
    if source.duplicated(keys).any():
        raise ValueError("Table 2 source has duplicate configured-group rows")
    counts = _group_count_frame(signed_rows)
    support = metadata[keys + ["support_size"]].copy()
    out = source[keys + ["rho_domain", "coherence_domain"]].merge(
        support, on=keys, how="left", validate="one_to_one"
    )
    out = out.merge(counts, on=keys, how="left", validate="one_to_one")
    if out.isna().any().any():
        raise ValueError("Table 2 join produced missing configured-group metadata")
    if not np.array_equal(out["selected_item_rows"], out["support_size"]):
        raise ValueError("Table 2 selected item-row counts differ from support sizes")
    out["configured_group"] = (
        out["domain"].astype(str) + " (" + out["chain"].astype(str) + ")"
    )
    return out[
        [
            "serotype",
            "configured_group",
            "support_size",
            "rho_domain",
            "coherence_domain",
            "selected_item_rows",
            "distinct_selected_fits_represented",
            "finite_signed_rows",
            "distinct_finite_fits_represented",
            "bh_flagged_rows",
            "distinct_bh_flagged_fits_represented",
        ]
    ].rename(
        columns={
            "serotype": "dataset",
            "rho_domain": "magnitude_rho",
            "coherence_domain": "operational_coherence",
            "bh_flagged_rows": "BH_flagged_rows",
            "distinct_bh_flagged_fits_represented": (
                "distinct_BH_flagged_fits_represented"
            ),
        }
    ).sort_values(["dataset", "configured_group"], kind="stable").reset_index(drop=True)


def _table3(paths: Paths, metadata: pd.DataFrame) -> pd.DataFrame:
    source = load_s7_figure(paths, "T3_catalytic_cross_serotype").copy()
    keys = ["serotype", "chain", "domain"]
    required = set(keys) | {"rho_domain", "beta_domain"}
    if missing := required - set(source.columns):
        raise ValueError(f"Table 3 source missing columns: {sorted(missing)}")
    support = metadata[keys + ["support_size"]]
    out = source.merge(support, on=keys, how="left", validate="one_to_one")
    if out["support_size"].isna().any():
        raise ValueError("Table 3 lacks configured-group support metadata")
    return out[
        ["serotype", "chain", "domain", "support_size", "rho_domain", "beta_domain"]
    ].rename(
        columns={
            "serotype": "dataset",
            "domain": "configured_group",
            "rho_domain": "magnitude_rho",
            "beta_domain": "fitted_magnitude_beta",
        }
    ).sort_values(
        ["configured_group", "dataset"], kind="stable"
    ).reset_index(drop=True)


def _table4_authoritative() -> tuple[pd.DataFrame, dict[str, int]]:
    summary = pd.read_csv(PHASE3E_DIR / "figure5_coordinate_summary.csv")
    half = summary[np.isclose(summary["cutoff"], PROVISIONAL_RHO, atol=1e-12)].copy()
    if len(half) != 248:
        raise ValueError(f"Table 4 expected 248 coordinates at 0.5, found {len(half)}")
    strict = half[half["crossing_all_four"].astype(bool)].copy()
    all_where_present = int((half["n_crossing"] == half["n_present"]).sum())
    if len(strict) != 18 or all_where_present != 30:
        raise ValueError(
            "Table 4 registered recurrence counts differ: "
            f"strict={len(strict)}, all-where-present={all_where_present}"
        )
    out = strict[
        [
            "canonical_coordinate",
            "chain",
            "configured_group_label",
            "n_total",
            "n_present",
            "n_absent",
            "n_crossing",
            "n_present_below_cutoff",
            "present_in_all_four",
            "crossing_all_four",
        ]
    ].rename(columns={"configured_group_label": "configured_group"})
    out["_coordinate_number"] = out["canonical_coordinate"].map(
        lambda value: parse_canon_label(value)[1]
    )
    out = (
        out.sort_values(
            ["chain", "_coordinate_number", "canonical_coordinate"], kind="stable"
        )
        .drop(columns="_coordinate_number")
        .reset_index(drop=True)
    )
    return out, {
        "union_coordinates": 248,
        "partial_presence_coordinates": int((half["n_present"] < 4).sum()),
        "all_where_present_count": all_where_present,
        "strict_all_four_count": len(strict),
    }


def _table4_generic(paths: Paths) -> tuple[pd.DataFrame, dict[str, int]]:
    position = load_table(paths, "s5", "position_conservation").copy()
    required = {
        "canon_label",
        "chain",
        "domain",
        "n_serotypes_present",
        "n_serotypes_reproducible",
    }
    if missing := required - set(position.columns):
        raise ValueError(f"generic Table 4 source missing columns: {sorted(missing)}")
    position["n_total"] = 4
    position["n_present"] = position["n_serotypes_present"].astype(int)
    position["n_absent"] = 4 - position["n_present"]
    position["n_crossing"] = position["n_serotypes_reproducible"].astype(int)
    position["n_present_below_cutoff"] = (
        position["n_present"] - position["n_crossing"]
    )
    position["present_in_all_four"] = position["n_present"] == 4
    position["crossing_all_four"] = (
        position["present_in_all_four"] & (position["n_crossing"] == 4)
    )
    strict = position[position["crossing_all_four"]].copy()
    out = strict[
        [
            "canon_label",
            "chain",
            "domain",
            "n_total",
            "n_present",
            "n_absent",
            "n_crossing",
            "n_present_below_cutoff",
            "present_in_all_four",
            "crossing_all_four",
        ]
    ].rename(
        columns={
            "canon_label": "canonical_coordinate",
            "domain": "configured_group",
        }
    )
    return out.reset_index(drop=True), {
        "union_coordinates": len(position),
        "partial_presence_coordinates": int((position["n_present"] < 4).sum()),
        "all_where_present_count": int(
            (position["n_crossing"] == position["n_present"]).sum()
        ),
        "strict_all_four_count": len(strict),
    }


def _table5(paths: Paths, metadata: pd.DataFrame) -> pd.DataFrame:
    source = load_s7_figure(paths, "T5_variance_component_budget").copy()
    keys = ["serotype", "chain", "domain"]
    required = set(keys) | {
        "tau2",
        "sigma2_bar",
        "frac_tau2",
        "frac_sigma2",
        "variance_regime",
    }
    if missing := required - set(source.columns):
        raise ValueError(f"Table 5 source missing columns: {sorted(missing)}")
    if source.duplicated(keys).any():
        raise ValueError("Table 5 source has duplicate configured-group rows")
    support = metadata[keys + ["support_size"]]
    out = source.merge(support, on=keys, how="left", validate="one_to_one")
    if out["support_size"].isna().any():
        raise ValueError("Table 5 lacks configured-group support metadata")
    out["configured_group"] = (
        out["domain"].astype(str) + " (" + out["chain"].astype(str) + ")"
    )
    out = out[
        [
            "serotype",
            "configured_group",
            "support_size",
            "tau2",
            "sigma2_bar",
            "frac_tau2",
            "frac_sigma2",
            "variance_regime",
        ]
    ].rename(
        columns={
            "serotype": "dataset",
            "tau2": "fitted_between_run_component",
            "sigma2_bar": "propagated_within_run_component",
            "frac_tau2": "between_run_fraction",
            "frac_sigma2": "within_run_fraction",
            "variance_regime": "operational_variance_regime",
        }
    ).sort_values(["dataset", "configured_group"], kind="stable").reset_index(drop=True)
    out["operational_variance_regime"] = out[
        "operational_variance_regime"
    ].astype(str).str.replace("_", "-", regex=False)
    return out


def prepare_tables(paths: Paths) -> PreparedTables:
    """Build all repaired table frames and validate cross-table invariants."""
    source_hashes = _source_hashes(paths)
    signed = load_signed_figure_data(paths)
    rows = signed.selected_rows
    metadata, metadata_provenance = _group_metadata(paths)
    table1 = _table1(paths, rows)
    table2_full = _table2(paths, rows, metadata)
    table2 = table2_full[
        [
            "dataset",
            "configured_group",
            "support_size",
            "magnitude_rho",
            "operational_coherence",
            "selected_item_rows",
            "distinct_selected_fits_represented",
        ]
    ].copy()
    table2_signed = table2_full[
        [
            "dataset",
            "configured_group",
            "finite_signed_rows",
            "distinct_finite_fits_represented",
            "BH_flagged_rows",
            "distinct_BH_flagged_fits_represented",
        ]
    ].copy()
    table3 = _table3(paths, metadata)
    if source_hashes["outputs_s7/F1_reproducibility_landscape.parquet"] == AUTHORITATIVE_F1_HASH:
        table4, recurrence = _table4_authoritative()
    else:
        table4, recurrence = _table4_generic(paths)
    table5 = _table5(paths, metadata)
    counts = bookkeeping_counts(signed)
    regime_counts = (
        table5["operational_variance_regime"].value_counts().sort_index().to_dict()
    )
    consistency = {
        **counts,
        **recurrence,
        "table5_rows": len(table5),
        "variance_regime_counts": {k: int(v) for k, v in regime_counts.items()},
        "metadata_provenance": metadata_provenance,
    }
    if signed.authoritative:
        expected = {
            "raw_interval_exclusion_rows": 181,
            "distinct_raw_interval_exclusion_fits": 174,
            "bh_flagged_rows": 152,
            "distinct_bh_flagged_fits": 145,
            "finite_signed_rows": 302,
            "distinct_finite_fits": 268,
        }
        for key, value in expected.items():
            if consistency[key] != value:
                raise ValueError(
                    f"cross-table invariant {key} expected {value}, "
                    f"observed {consistency[key]}"
                )
        if recurrence != {
            "union_coordinates": 248,
            "partial_presence_coordinates": 49,
            "all_where_present_count": 30,
            "strict_all_four_count": 18,
        }:
            raise ValueError(f"Table 4 recurrence invariants differ: {recurrence}")
        if consistency["variance_regime_counts"] != {
            "balanced": 2,
            "replicate-dominated": 28,
            "sampling-dominated": 2,
        }:
            raise ValueError("Table 5 28/2/2 classification differs")
    return PreparedTables(
        frames={
            "table1_per_serotype_summary": table1,
            "table2_configured_group_summary": table2,
            "supplementary_table2_conditional_signed_bookkeeping": table2_signed,
            "table3_configured_catalytic_groups": table3,
            "table4_all_four_operational_recurrence": table4,
            "table5_configured_group_variance_composition": table5,
        },
        source_hashes=source_hashes,
        consistency=consistency,
    )


def _assert_safe_table_text() -> None:
    text = " ".join(
        item
        for spec in TABLE_SPECS.values()
        for item in (spec["caption"], spec["note"])
    ).lower()
    for term in UNSAFE_TABLE_TERMS:
        if term in text:
            raise ValueError(f"unsafe table-facing term present: {term}")


def build(paths: Paths) -> list[Path]:
    """Export repaired CSV/LaTeX/HTML/Markdown tables and deterministic audits."""
    _assert_safe_table_text()
    prepared = prepare_tables(paths)
    written: list[Path] = []
    for stem, frame in prepared.frames.items():
        spec = TABLE_SPECS[stem]
        csv_path = paths.output / f"{stem}.csv"
        frame.to_csv(csv_path, index=False, lineterminator="\n")
        written.append(csv_path)
        written += export_table(
            frame,
            paths,
            stem,
            caption=spec["caption"],
            label=spec["label"],
            note=spec["note"],
        )

    output_hashes = {
        path.name: {"rows": len(prepared.frames[path.stem]), "sha256": sha256(path)}
        for path in written
        if path.suffix == ".csv"
    }
    manifest = {
        "schema_version": "phase3k-v1",
        "purpose": "presentation-only descriptive table repair",
        "script": {
            "path": str(Path(__file__).resolve()),
            "sha256": sha256(Path(__file__).resolve()),
        },
        "source_hashes": prepared.source_hashes,
        "output_csvs": output_hashes,
        "display_precision": "three significant digits in Markdown/HTML/LaTeX",
        "csv_precision": "pandas deterministic full stored representation",
        "provisional_rho_reference": PROVISIONAL_RHO,
        "operational_coherence_reference": OPERATIONAL_COHERENCE,
    }
    write_json(paths.output / "phase3k_table_manifest.json", manifest)
    write_json(
        paths.output / "phase3k_cross_table_consistency_audit.json",
        prepared.consistency,
    )
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for path in build(resolve_paths()):
        print(path)
