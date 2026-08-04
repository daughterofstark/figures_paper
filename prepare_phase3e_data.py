"""Prepare deterministic presentation-only data tables for Phase 3E.

Inputs
------
Frozen STRIDE downstream tables:

* ``outputs/stride_table.parquet``
* ``outputs/replicate_table.parquet``
* ``outputs_s2/residue_landscape.parquet``
* ``outputs_s5/position_conservation.parquet``

Completed V9A and V9B validation artifacts are read without modification.

Outputs
-------
CSV tables and a deterministic provenance manifest under
``derived_data/phase3e`` by default. This module does not plot, refit an
estimator, alter a threshold, or overwrite an authoritative input.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Iterable, Mapping, Sequence

import numpy as np
import pandas as pd

from figure_config import resolve_external_sources


STRIDE_COMMIT = "b06f7d1d40b5d860915cf1e1f046a6392f71a09b"
DOWNSTREAM_COMMIT = "4388dbcefaaa1672b63cc7bbce357e1b4ec40b71"
PRESENTATION_BASE_COMMIT = "b617444d9a3b5ef528a9d0cd1c9fc3dd738daadb"

DATASETS = ("DENV1", "DENV2", "DENV3", "DENV4")
DATASET_ID_MAP = {
    "dataset_1": "DENV1",
    "dataset_2": "DENV2",
    "dataset_3": "DENV3",
    "dataset_4": "DENV4",
}
NAMED_LEVELS = (
    "residue",
    "secondary_structure",
    "motif",
    "domain",
    "chain",
    "protein",
    "complex",
)
PRESENTATION_LEVELS = (
    "item",
    "configured_group",
    "distinct_chain",
    "root",
)
CROSSING_LEVELS = PRESENTATION_LEVELS + ("unresolved",)
CUTOFFS = (0.4, 0.5, 0.5330838393, 0.6)
TOP_N_VALUES = (5, 10, 20)
CATALYTIC_GROUPS = frozenset({"Catalytic Triad", "Oxyanion Loop"})
EXPECTED_DATASET_ITEMS = {
    "DENV1": 236,
    "DENV2": 247,
    "DENV3": 242,
    "DENV4": 199,
}
EXPECTED_FIGURE5_VECTORS = {
    0.4: (5, 26, 73, 91, 53),
    0.5: (18, 49, 100, 63, 18),
    0.5330838393: (26, 58, 98, 53, 13),
    0.6: (42, 84, 82, 33, 7),
}

OUTPUT_SORT_KEYS: Mapping[str, Sequence[str]] = {
    "figure2_unique_supports.csv": (
        "dataset",
        "canonical_item_key",
        "presentation_level_index",
    ),
    "figure2_item_first_aggregate.csv": ("dataset", "canonical_item_key"),
    "figure2_observed_first_crossing_sensitivity.csv": (
        "cutoff",
        "dataset",
        "crossing_level_index",
    ),
    "figure2_v9a_support_size.csv": (
        "dataset",
        "regime",
        "support_size_stratum_index",
    ),
    "figure2_v9a_support_size_association.csv": ("dataset", "regime"),
    "figure2_v9a_item_first_aggregate.csv": ("dataset", "regime"),
    "figure2_v9a_first_crossing.csv": (
        "dataset",
        "regime",
        "crossing_level_index",
    ),
    "figure2_v9a_alias_provenance.csv": ("metric"),
    "figure4_top1.csv": ("dataset", "run_index"),
    "figure4_topn_jaccard.csv": (
        "dataset",
        "N",
        "run_a_index",
        "run_b_index",
    ),
    "figure4_top10_heatmap.csv": ("dataset", "union_order", "run_index"),
    "figure5_status_matrix.csv": (
        "cutoff",
        "canonical_coordinate",
        "dataset",
    ),
    "figure5_coordinate_summary.csv": ("cutoff", "canonical_coordinate"),
    "figure5_row_order.csv": ("row_index",),
    "configured_group_supports.csv": (
        "dataset",
        "chain",
        "configured_group_label",
    ),
}


@dataclass(frozen=True)
class SourcePaths:
    """Authoritative input locations."""

    frozen_root: Path
    v9a_dir: Path
    v9b_dir: Path

    @property
    def stride_table(self) -> Path:
        return self.frozen_root / "outputs" / "stride_table.parquet"

    @property
    def replicate_table(self) -> Path:
        return self.frozen_root / "outputs" / "replicate_table.parquet"

    @property
    def residue_landscape(self) -> Path:
        return self.frozen_root / "outputs_s2" / "residue_landscape.parquet"

    @property
    def position_conservation(self) -> Path:
        return self.frozen_root / "outputs_s5" / "position_conservation.parquet"

    @property
    def v9a_support_size(self) -> Path:
        return self.v9a_dir / "v9a_frozen_design_support_size.csv"

    @property
    def v9a_first_crossing(self) -> Path:
        return self.v9a_dir / "v9a_frozen_design_first_crossing.csv"

    @property
    def v9a_alias_effects(self) -> Path:
        return self.v9a_dir / "v9a_frozen_design_alias_effects.csv"

    @property
    def v9a_results(self) -> Path:
        return self.v9a_dir / "v9a_frozen_design_results.json"

    @property
    def v9a_manifest(self) -> Path:
        return self.v9a_dir / "v9a_frozen_design_manifest.json"

    @property
    def v9a_snapshot(self) -> Path:
        return self.v9a_dir / "v9a_frozen_design_input.json"

    @property
    def v9b_results(self) -> Path:
        return self.v9b_dir / "v9b_partial_presence_results.json"

    @property
    def v9b_fixture(self) -> Path:
        return self.v9b_dir / "v9b_partial_presence_fixture.csv"

    @property
    def v9b_expected(self) -> Path:
        return self.v9b_dir / "v9b_partial_presence_expected.csv"

    @property
    def v9b_manifest(self) -> Path:
        return self.v9b_dir / "v9b_manifest.json"

    def all_inputs(self) -> tuple[Path, ...]:
        return (
            self.stride_table,
            self.replicate_table,
            self.residue_landscape,
            self.position_conservation,
            self.v9a_support_size,
            self.v9a_first_crossing,
            self.v9a_alias_effects,
            self.v9a_results,
            self.v9a_manifest,
            self.v9a_snapshot,
            self.v9b_results,
            self.v9b_fixture,
            self.v9b_expected,
            self.v9b_manifest,
        )


def sha256_file(path: Path) -> str:
    """Return a file's SHA-256 digest."""

    h = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def membership_hash(members: Iterable[str]) -> str:
    """Hash a sorted exact support membership."""

    payload = "\n".join(sorted(str(v) for v in members)).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def canonical_sort_key(value: str) -> tuple[str, int, str]:
    """Sort ``CHAIN:number`` keys deterministically."""

    text = str(value)
    if ":" not in text:
        return (text, 0, text)
    chain, residue = text.rsplit(":", 1)
    try:
        number = int(residue)
    except ValueError:
        number = 0
    return (chain, number, text)


def _require_columns(df: pd.DataFrame, required: Sequence[str], source: Path) -> None:
    missing = sorted(set(required) - set(df.columns))
    if missing:
        raise ValueError(
            f"{source} is missing required columns {missing}; "
            f"available columns are {list(df.columns)}"
        )


def _load_json(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        value = json.load(handle)
    if not isinstance(value, dict):
        raise ValueError(f"expected a JSON object in {path}")
    return value


def _represent_members(members: Iterable[str]) -> str:
    return ";".join(sorted((str(v) for v in members), key=canonical_sort_key))


def _base_item_annotations(replicates: pd.DataFrame, source: Path) -> pd.DataFrame:
    required = (
        "serotype",
        "replicate_index",
        "canon_label",
        "chain",
        "domain",
        "region_id",
        "r",
    )
    _require_columns(replicates, required, source)
    keys = ["serotype", "canon_label"]
    annotation_columns = ["chain", "domain", "region_id"]
    for column in annotation_columns:
        counts = replicates.groupby(keys, sort=True)[column].nunique(dropna=False)
        if (counts != 1).any():
            bad = counts[counts != 1].index.tolist()[:5]
            raise ValueError(f"non-constant {column} across runs for {bad}")
    items = (
        replicates[keys + annotation_columns]
        .drop_duplicates(keys)
        .rename(columns={"serotype": "dataset", "canon_label": "canonical_item_key"})
        .sort_values(["dataset", "canonical_item_key"], kind="stable")
        .reset_index(drop=True)
    )
    observed = items.groupby("dataset").size().to_dict()
    if observed != EXPECTED_DATASET_ITEMS:
        raise ValueError(f"unexpected dataset item counts: {observed}")
    return items


def _support_lookup(
    items: pd.DataFrame,
    group_column: str,
) -> dict[tuple[str, str], tuple[str, ...]]:
    lookup: dict[tuple[str, str], tuple[str, ...]] = {}
    for (dataset, label), group in items.groupby(
        ["dataset", group_column], sort=True, dropna=False
    ):
        members = tuple(
            sorted(group["canonical_item_key"].astype(str), key=canonical_sort_key)
        )
        lookup[(str(dataset), str(label))] = members
    return lookup


def prepare_unique_supports(
    stride: pd.DataFrame,
    replicates: pd.DataFrame,
    *,
    stride_source: Path,
    replicate_source: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Collapse exact hierarchy aliases into unique presentation supports."""

    _require_columns(
        stride,
        (
            "serotype",
            "canon_label",
            "scale_level",
            "scale_index",
            "rho",
            "h_chain",
            "h_domain",
        ),
        stride_source,
    )
    items = _base_item_annotations(replicates, replicate_source)
    chain_members = _support_lookup(items, "chain")
    group_members = _support_lookup(items, "region_id")
    root_members = {
        str(dataset): tuple(
            sorted(group["canonical_item_key"].astype(str), key=canonical_sort_key)
        )
        for dataset, group in items.groupby("dataset", sort=True)
    }
    named_order = {level: index for index, level in enumerate(NAMED_LEVELS)}
    presentation_order = {
        level: index for index, level in enumerate(PRESENTATION_LEVELS)
    }
    stride = stride.copy()
    stride["serotype"] = stride["serotype"].astype(str)
    stride["canon_label"] = stride["canon_label"].astype(str)
    stride["scale_level"] = stride["scale_level"].astype(str)
    stride["rho"] = pd.to_numeric(stride["rho"], errors="coerce")
    if not np.isfinite(stride["rho"]).all():
        raise ValueError(f"non-finite rho in {stride_source}")
    profile_counts = stride.groupby(["serotype", "canon_label"]).size()
    if len(profile_counts) != 924 or set(profile_counts.unique()) != {7}:
        raise ValueError(
            "stride table must contain 924 profiles with seven named rows each"
        )

    item_meta = items.set_index(["dataset", "canonical_item_key"])
    rows: list[dict[str, object]] = []
    ss_domain_aliases = 0
    root_aliases = 0
    group_chain_aliases = 0

    for (dataset, item), profile in stride.groupby(
        ["serotype", "canon_label"], sort=True
    ):
        profile = profile.sort_values("scale_index", kind="stable")
        if tuple(profile["scale_level"]) != NAMED_LEVELS:
            raise ValueError(f"unexpected named hierarchy for {(dataset, item)}")
        meta = item_meta.loc[(dataset, item)]
        named_supports: list[tuple[str, tuple[str, ...], float]] = []
        for row in profile.itertuples(index=False):
            level = str(row.scale_level)
            if level == "residue":
                members = (str(item),)
            elif level in {"secondary_structure", "motif", "domain"}:
                members = group_members[(str(dataset), str(meta.region_id))]
            elif level == "chain":
                members = chain_members[(str(dataset), str(meta.chain))]
            elif level in {"protein", "complex"}:
                members = root_members[str(dataset)]
            else:  # pragma: no cover - guarded by named hierarchy check
                raise ValueError(f"unknown hierarchy level {level}")
            named_supports.append((level, members, float(row.rho)))

        grouped: dict[tuple[str, ...], list[tuple[str, float]]] = {}
        for level, members, rho in named_supports:
            grouped.setdefault(members, []).append((level, rho))

        if [
            level for level, _ in grouped[
                group_members[(str(dataset), str(meta.region_id))]
            ]
            if level in {"secondary_structure", "motif", "domain"}
        ] != ["secondary_structure", "motif", "domain"]:
            raise ValueError(f"configured-group aliases differ for {(dataset, item)}")
        ss_domain_aliases += 1
        root_levels = [
            level
            for level, _ in grouped[root_members[str(dataset)]]
            if level in {"protein", "complex"}
        ]
        if root_levels != ["protein", "complex"]:
            raise ValueError(f"root aliases differ for {(dataset, item)}")
        root_aliases += 1

        collapsed_rows: list[dict[str, object]] = []
        for members, representations in grouped.items():
            representations = sorted(
                representations, key=lambda value: named_order[value[0]]
            )
            levels = [level for level, _ in representations]
            rho_values = np.asarray([rho for _, rho in representations], dtype=float)
            if np.max(np.abs(rho_values - rho_values[0])) > 1e-12:
                raise ValueError(
                    f"rho differs across exact aliases for {(dataset, item, levels)}"
                )
            if "residue" in levels:
                presentation_level = "item"
            elif any(
                level in {"secondary_structure", "motif", "domain"}
                for level in levels
            ):
                presentation_level = "configured_group"
            elif "chain" in levels:
                presentation_level = "distinct_chain"
            elif any(level in {"protein", "complex"} for level in levels):
                presentation_level = "root"
            else:  # pragma: no cover
                raise ValueError(f"cannot classify support levels {levels}")
            if (
                presentation_level == "configured_group"
                and "chain" in levels
            ):
                group_chain_aliases += 1
            member_text = _represent_members(members)
            digest = membership_hash(members)
            collapsed_rows.append(
                {
                    "dataset": str(dataset),
                    "canonical_item_key": str(item),
                    "chain": str(meta.chain),
                    "configured_group_label": str(meta.domain),
                    "original_named_levels": "|".join(levels),
                    "collapsed_presentation_level": presentation_level,
                    "presentation_level_index": presentation_order[presentation_level],
                    "exact_support_identifier": f"{dataset}:{digest}",
                    "support_size": len(members),
                    "rho": float(rho_values[0]),
                    "support_membership_hash": digest,
                    "support_membership": member_text,
                    "is_singleton": presentation_level == "item",
                    "is_first_unique_aggregate": False,
                    "is_distinct_chain_expansion": presentation_level
                    == "distinct_chain",
                    "is_whole_system_root": presentation_level == "root",
                }
            )
        collapsed_rows.sort(key=lambda row: int(row["presentation_level_index"]))
        aggregate_indices = [
            index for index, row in enumerate(collapsed_rows) if not row["is_singleton"]
        ]
        if not aggregate_indices:
            raise ValueError(f"no unique aggregate for {(dataset, item)}")
        collapsed_rows[aggregate_indices[0]]["is_first_unique_aggregate"] = True
        rows.extend(collapsed_rows)

    unique = pd.DataFrame(rows).sort_values(
        ["dataset", "canonical_item_key", "presentation_level_index"],
        kind="stable",
    )
    singletons = unique[unique["is_singleton"]]
    first = unique[unique["is_first_unique_aggregate"]]
    if len(singletons) != 924 or len(first) != 924:
        raise ValueError(
            f"expected 924 singleton and aggregate rows, got "
            f"{len(singletons)} and {len(first)}"
        )
    paired = singletons[
        ["dataset", "canonical_item_key", "rho"]
    ].rename(columns={"rho": "item_rho"}).merge(
        first[
            [
                "dataset",
                "canonical_item_key",
                "rho",
                "support_size",
                "support_membership_hash",
            ]
        ].rename(
            columns={
                "rho": "first_aggregate_rho",
                "support_size": "first_aggregate_support_size",
                "support_membership_hash": "first_aggregate_membership_hash",
            }
        ),
        on=["dataset", "canonical_item_key"],
        validate="one_to_one",
    )
    paired["delta_rho"] = paired["first_aggregate_rho"] - paired["item_rho"]
    paired = paired.sort_values(
        ["dataset", "canonical_item_key"], kind="stable"
    ).reset_index(drop=True)
    diagnostics = {
        "item_profiles": int(len(singletons)),
        "unique_item_support_rows": int(len(unique)),
        "ss_motif_domain_alias_profiles": int(ss_domain_aliases),
        "protein_complex_alias_profiles": int(root_aliases),
        "group_chain_alias_profiles": int(group_chain_aliases),
    }
    expected_diagnostics = {
        "item_profiles": 924,
        "ss_motif_domain_alias_profiles": 924,
        "protein_complex_alias_profiles": 924,
        "group_chain_alias_profiles": 203,
    }
    for key, expected in expected_diagnostics.items():
        if diagnostics[key] != expected:
            raise ValueError(f"{key}: expected {expected}, got {diagnostics[key]}")
    return unique.reset_index(drop=True), paired, diagnostics


def prepare_observed_crossing_sensitivity(unique: pd.DataFrame) -> pd.DataFrame:
    """Count first unique-support crossings at predeclared cutoffs."""

    rows: list[dict[str, object]] = []
    datasets = list(DATASETS)
    for cutoff in CUTOFFS:
        per_item: list[tuple[str, str]] = []
        for (dataset, item), profile in unique.groupby(
            ["dataset", "canonical_item_key"], sort=True
        ):
            profile = profile.sort_values("presentation_level_index", kind="stable")
            crossed = profile[profile["rho"] >= cutoff]
            level = (
                str(crossed.iloc[0]["collapsed_presentation_level"])
                if not crossed.empty
                else "unresolved"
            )
            per_item.append((str(dataset), level))
        decisions = pd.DataFrame(per_item, columns=["dataset", "crossing_level"])
        for dataset in datasets + ["ALL"]:
            subset = decisions if dataset == "ALL" else decisions[decisions["dataset"] == dataset]
            denominator = len(subset)
            for index, level in enumerate(CROSSING_LEVELS):
                count = int((subset["crossing_level"] == level).sum())
                rows.append(
                    {
                        "cutoff": cutoff,
                        "dataset": dataset,
                        "crossing_level": level,
                        "crossing_level_index": index,
                        "count": count,
                        "denominator": denominator,
                        "fraction": count / denominator,
                    }
                )
    return pd.DataFrame(rows).sort_values(
        ["cutoff", "dataset", "crossing_level_index"], kind="stable"
    )


def _map_v9a_dataset(frame: pd.DataFrame) -> pd.DataFrame:
    frame = frame.copy()
    frame["dataset"] = frame["dataset_id"].map(DATASET_ID_MAP)
    if frame["dataset"].isna().any():
        raise ValueError("unknown V9A dataset_id")
    return frame


def prepare_v9a_tables(
    paths: SourcePaths,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Extract presentation summaries from completed V9A artifacts."""

    v9a_manifest = _load_json(paths.v9a_manifest)
    recorded = v9a_manifest.get("artifact_sha256")
    if not isinstance(recorded, dict):
        raise ValueError("V9A manifest lacks artifact_sha256")
    expected_artifacts = {
        "snapshot": paths.v9a_snapshot,
        "results": paths.v9a_results,
        "support_size": paths.v9a_support_size,
        "first_crossing": paths.v9a_first_crossing,
        "alias_effects": paths.v9a_alias_effects,
    }
    for key, path in expected_artifacts.items():
        actual = sha256_file(path)
        if recorded.get(key) != actual:
            raise ValueError(
                f"V9A hash mismatch for {key}: {recorded.get(key)} != {actual}"
            )

    support = pd.read_csv(paths.v9a_support_size)
    support = support[support["group_type"] == "stratum"].copy()
    support = _map_v9a_dataset(support)
    stratum_order = {
        value: index
        for index, value in enumerate(
            ("1", "2-4", "5-8", "9-16", "17-32", "33-64", "65-128", "129+")
        )
    }
    support["support_size_stratum"] = support["support_size_group"].astype(str)
    support["support_size_stratum_index"] = support[
        "support_size_stratum"
    ].map(stratum_order)
    if support["support_size_stratum_index"].isna().any():
        raise ValueError("unknown V9A support-size stratum")
    support_columns = [
        "dataset",
        "regime",
        "support_size_stratum",
        "support_size_stratum_index",
        "distinct_supports",
        "expected_fit_count",
        "finite_fit_count",
        "mean_rho",
        "median_rho",
        "across_draw_mean_rho",
        "mc_t95_lower",
        "mc_t95_upper",
        "distribution_q025",
        "distribution_q25",
        "distribution_q50",
        "distribution_q75",
        "distribution_q975",
    ]
    support = support[support_columns].sort_values(
        ["dataset", "regime", "support_size_stratum_index"], kind="stable"
    )

    results = _load_json(paths.v9a_results)
    association_value = results.get("support_size_association")
    if not isinstance(association_value, list):
        raise ValueError("V9A results lack support_size_association")
    association = _map_v9a_dataset(pd.DataFrame(association_value))
    association = association[
        [
            "dataset",
            "regime",
            "expected_coefficients",
            "n_coefficients",
            "supports_per_coefficient_min",
            "supports_per_coefficient_max",
            "mean_coefficient",
            "mc_t95_lower",
            "mc_t95_upper",
            "median",
            "distribution_q025",
            "distribution_q975",
        ]
    ].rename(columns={"mean_coefficient": "mean_spearman"})
    association = association.sort_values(["dataset", "regime"], kind="stable")

    aggregate_value = results.get("item_to_first_unique_aggregate")
    if not isinstance(aggregate_value, list):
        raise ValueError("V9A results lack item_to_first_unique_aggregate")
    aggregate_rows: list[dict[str, object]] = []
    for record in aggregate_value:
        if not isinstance(record, dict) or not isinstance(record.get("summary"), dict):
            raise ValueError("malformed V9A item-to-aggregate record")
        summary = record["summary"]
        row: dict[str, object] = {
            "dataset_id": record["dataset_id"],
            "regime": record["regime"],
            "expected_draws": record["expected_draws"],
            "successful_draws": record["successful_draws"],
            "failed_draws": record["failed_draws"],
            "item_comparisons_per_draw": record["item_comparisons_per_draw"],
            "total_expected_comparisons": record["total_expected_comparisons"],
        }
        for metric in (
            "fraction_increased",
            "fraction_equal",
            "fraction_decreased",
            "median_delta_rho",
        ):
            metric_value = summary.get(metric)
            if not isinstance(metric_value, dict):
                raise ValueError(f"V9A aggregate summary lacks {metric}")
            row[f"{metric}_mean"] = metric_value["mean"]
            row[f"{metric}_mc_t95_lower"] = metric_value["mc_t95_lower"]
            row[f"{metric}_mc_t95_upper"] = metric_value["mc_t95_upper"]
        aggregate_rows.append(row)
    aggregate = _map_v9a_dataset(pd.DataFrame(aggregate_rows))
    aggregate = aggregate.drop(columns="dataset_id").sort_values(
        ["dataset", "regime"], kind="stable"
    )

    crossing = pd.read_csv(paths.v9a_first_crossing)
    crossing = crossing[
        crossing["hierarchy_variant"] == "collapsed_unique_support"
    ].copy()
    crossing = crossing[
        crossing["crossing_level"].isin(
            ["item", "configured_group", "chain", "whole_system_root", "unresolved"]
        )
    ]
    crossing["crossing_level"] = crossing["crossing_level"].replace(
        {"chain": "distinct_chain", "whole_system_root": "root"}
    )
    crossing["crossing_level_index"] = crossing["crossing_level"].map(
        {level: index for index, level in enumerate(CROSSING_LEVELS)}
    )
    crossing = _map_v9a_dataset(crossing)
    crossing = crossing[
        [
            "dataset",
            "regime",
            "crossing_level",
            "crossing_level_index",
            "expected_draws",
            "successful_draws",
            "failed_draws",
            "item_profiles_per_draw",
            "expected_item_decisions",
            "total_crossing_count_successful_draws",
            "mean_draw_count",
            "mean_draw_fraction",
            "mc_t95_lower",
            "mc_t95_upper",
        ]
    ].sort_values(
        ["dataset", "regime", "crossing_level_index"], kind="stable"
    )

    alias_source = pd.read_csv(paths.v9a_alias_effects)
    paired = alias_source[alias_source["record_type"] == "paired_metric"]
    total_decisions = int(
        paired.groupby(["dataset_id", "regime"])["successful_item_decisions"]
        .first()
        .sum()
    )
    metric_map = {
        "reported_label_changed": "labels_changed",
        "selected_membership_changed": "selected_membership_changed",
        "selected_rho_changed_gt_1e-12": "selected_rho_changed_gt_1e-12",
    }
    alias_rows = [
        {
            "metric": "total_decisions",
            "count": total_decisions,
            "fraction": 1.0,
        }
    ]
    for source_metric, display_metric in metric_map.items():
        count = int(
            paired.loc[
                paired["metric_or_alias_block"] == source_metric, "count"
            ].sum()
        )
        alias_rows.append(
            {
                "metric": display_metric,
                "count": count,
                "fraction": count / total_decisions,
            }
        )
    alias = pd.DataFrame(alias_rows)
    alias["_order"] = alias["metric"].map(
        {
            "total_decisions": 0,
            "labels_changed": 1,
            "selected_membership_changed": 2,
            "selected_rho_changed_gt_1e-12": 3,
        }
    )
    alias = alias.sort_values("_order", kind="stable").drop(columns="_order")
    expected_alias = {
        "total_decisions": 924000,
        "labels_changed": 672379,
        "selected_membership_changed": 0,
        "selected_rho_changed_gt_1e-12": 0,
    }
    observed_alias = dict(zip(alias["metric"], alias["count"], strict=True))
    if observed_alias != expected_alias:
        raise ValueError(f"unexpected V9A alias summary: {observed_alias}")
    return (
        support.reset_index(drop=True),
        association.reset_index(drop=True),
        aggregate.reset_index(drop=True),
        crossing.reset_index(drop=True),
        alias.reset_index(drop=True),
    )


def _status_for(rho: float | None, cutoff: float) -> str:
    if rho is None or not np.isfinite(rho):
        return "ABSENT"
    return "CROSSING" if rho >= cutoff else "PRESENT_BELOW"


def prepare_figure5_tables(
    residue: pd.DataFrame,
    positions: pd.DataFrame,
    *,
    residue_source: Path,
    position_source: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Prepare presence-aware status, summary, and row-order tables."""

    _require_columns(
        residue,
        ("serotype", "canon_label", "rho_residue", "chain", "domain"),
        residue_source,
    )
    _require_columns(
        positions,
        ("canon_label", "chain", "domain", "n_serotypes_present"),
        position_source,
    )
    residue = residue.copy()
    residue["rho_residue"] = pd.to_numeric(residue["rho_residue"], errors="coerce")
    if residue.duplicated(["serotype", "canon_label"]).any():
        raise ValueError("duplicate item rows in residue landscape")
    union = (
        positions[["canon_label", "chain", "domain"]]
        .drop_duplicates("canon_label")
        .sort_values("canon_label", kind="stable")
    )
    if len(union) != 248:
        raise ValueError(f"expected 248 union coordinates, got {len(union)}")
    rho_lookup = {
        (str(row.serotype), str(row.canon_label)): float(row.rho_residue)
        for row in residue.itertuples(index=False)
    }
    status_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []
    for cutoff in CUTOFFS:
        for coordinate, chain, group_label in union.itertuples(index=False, name=None):
            states: list[str] = []
            for dataset in DATASETS:
                rho = rho_lookup.get((dataset, str(coordinate)))
                state = _status_for(rho, cutoff)
                states.append(state)
                status_rows.append(
                    {
                        "cutoff": cutoff,
                        "canonical_coordinate": str(coordinate),
                        "dataset": dataset,
                        "chain": str(chain),
                        "configured_group_label": str(group_label),
                        "rho": rho if rho is not None else np.nan,
                        "state": state,
                    }
                )
            n_present = sum(state != "ABSENT" for state in states)
            n_crossing = states.count("CROSSING")
            summary_rows.append(
                {
                    "cutoff": cutoff,
                    "canonical_coordinate": str(coordinate),
                    "chain": str(chain),
                    "configured_group_label": str(group_label),
                    "n_total": 4,
                    "n_present": n_present,
                    "n_absent": 4 - n_present,
                    "n_crossing": n_crossing,
                    "n_present_below_cutoff": n_present - n_crossing,
                    "present_in_all_four": n_present == 4,
                    "crossing_all_four": n_present == 4 and n_crossing == 4,
                }
            )
    status = pd.DataFrame(status_rows).sort_values(
        ["cutoff", "canonical_coordinate", "dataset"], kind="stable"
    )
    summary = pd.DataFrame(summary_rows).sort_values(
        ["cutoff", "canonical_coordinate"], kind="stable"
    )
    at_half = summary[summary["cutoff"] == 0.5]
    partial_presence = int((at_half["n_present"] < 4).sum())
    if partial_presence != 49:
        raise ValueError(
            f"expected 49 partial-presence coordinates, got {partial_presence}"
        )
    for cutoff, expected in EXPECTED_FIGURE5_VECTORS.items():
        subset = summary[summary["cutoff"] == cutoff]
        vector = tuple(
            int(value)
            for value in subset["n_crossing"]
            .value_counts()
            .reindex(range(5), fill_value=0)
            .to_numpy()
        )
        if vector != expected:
            raise ValueError(
                f"Figure 5 vector mismatch at {cutoff}: {vector} != {expected}"
            )

    order_base = at_half.copy()
    order_base["_chain_sort"] = order_base["chain"].astype(str)
    order_base["_coordinate_number"] = order_base["canonical_coordinate"].map(
        lambda value: canonical_sort_key(value)[1]
    )
    order_base = order_base.sort_values(
        [
            "n_present",
            "n_crossing",
            "_chain_sort",
            "_coordinate_number",
            "canonical_coordinate",
        ],
        ascending=[False, False, True, True, True],
        kind="stable",
    ).reset_index(drop=True)
    order_base.insert(0, "row_index", np.arange(len(order_base), dtype=int))
    order = order_base[
        [
            "row_index",
            "canonical_coordinate",
            "chain",
            "configured_group_label",
            "n_present",
            "n_crossing",
        ]
    ]
    return status.reset_index(drop=True), summary.reset_index(drop=True), order


def _rank_effects(replicates: pd.DataFrame, source: Path) -> tuple[pd.DataFrame, dict[str, int]]:
    required = (
        "serotype",
        "replicate",
        "replicate_index",
        "canon_label",
        "r",
        "chain",
        "domain",
    )
    _require_columns(replicates, required, source)
    frame = replicates[list(required)].copy()
    frame = frame.rename(
        columns={
            "serotype": "dataset",
            "replicate": "run",
            "replicate_index": "run_index",
            "canon_label": "canonical_coordinate",
            "r": "effect",
            "domain": "configured_group_label",
        }
    )
    frame["effect"] = pd.to_numeric(frame["effect"], errors="coerce")
    if not np.isfinite(frame["effect"]).all():
        raise ValueError("ranking baseline contains non-finite effects")
    frame["absolute_effect"] = frame["effect"].abs()
    frame["sign"] = np.sign(frame["effect"]).astype(int)
    ranked_groups: list[pd.DataFrame] = []
    tie_groups = 0
    tied_rows = 0
    for _, group in frame.groupby(["dataset", "run_index"], sort=True):
        counts = group["absolute_effect"].value_counts()
        tie_groups += int((counts > 1).sum())
        tied_rows += int(counts[counts > 1].sum())
        ordered = group.sort_values(
            ["absolute_effect", "canonical_coordinate"],
            ascending=[False, True],
            kind="stable",
        ).copy()
        ordered["rank"] = np.arange(1, len(ordered) + 1, dtype=int)
        ranked_groups.append(ordered)
    ranked = pd.concat(ranked_groups, ignore_index=True)
    return ranked, {"exact_tie_groups": tie_groups, "rows_in_exact_ties": tied_rows}


def prepare_figure4_tables(
    replicates: pd.DataFrame,
    *,
    replicate_source: Path,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict[str, int]]:
    """Prepare deterministic top-1, Jaccard, and top-10 union tables."""

    ranked, tie_diagnostics = _rank_effects(replicates, replicate_source)
    top1 = ranked[ranked["rank"] == 1][
        [
            "dataset",
            "run",
            "run_index",
            "canonical_coordinate",
            "chain",
            "configured_group_label",
            "effect",
            "sign",
            "absolute_effect",
        ]
    ].rename(columns={"canonical_coordinate": "top_item"})
    top1["canonical_coordinate"] = top1["top_item"]
    top1 = top1[
        [
            "dataset",
            "run",
            "run_index",
            "top_item",
            "canonical_coordinate",
            "chain",
            "configured_group_label",
            "effect",
            "sign",
            "absolute_effect",
        ]
    ].sort_values(["dataset", "run_index"], kind="stable")
    if len(top1) != 12:
        raise ValueError(f"expected 12 top-1 rows, got {len(top1)}")

    jaccard_rows: list[dict[str, object]] = []
    for dataset, dataset_frame in ranked.groupby("dataset", sort=True):
        runs = (
            dataset_frame[["run", "run_index"]]
            .drop_duplicates()
            .sort_values("run_index")
        )
        run_pairs = list(runs.itertuples(index=False, name=None))
        if len(run_pairs) != 3:
            raise ValueError(f"{dataset} does not have exactly three runs")
        for n in TOP_N_VALUES:
            sets = {
                int(run_index): set(
                    dataset_frame[
                        (dataset_frame["run_index"] == run_index)
                        & (dataset_frame["rank"] <= n)
                    ]["canonical_coordinate"]
                )
                for _, run_index in run_pairs
            }
            labels = {int(index): str(run) for run, index in run_pairs}
            for run_a_index, run_b_index in combinations(sorted(sets), 2):
                set_a = sets[run_a_index]
                set_b = sets[run_b_index]
                union = set_a | set_b
                jaccard_rows.append(
                    {
                        "dataset": str(dataset),
                        "N": n,
                        "run_a": labels[run_a_index],
                        "run_a_index": run_a_index,
                        "run_b": labels[run_b_index],
                        "run_b_index": run_b_index,
                        "intersection_count": len(set_a & set_b),
                        "union_count": len(union),
                        "jaccard": len(set_a & set_b) / len(union),
                    }
                )
    jaccard = pd.DataFrame(jaccard_rows).sort_values(
        ["dataset", "N", "run_a_index", "run_b_index"], kind="stable"
    )
    if len(jaccard) != 36:
        raise ValueError(f"expected 36 pairwise top-N rows, got {len(jaccard)}")

    heatmap_rows: list[dict[str, object]] = []
    for dataset, dataset_frame in ranked.groupby("dataset", sort=True):
        top10 = dataset_frame[dataset_frame["rank"] <= 10]
        union_items = sorted(
            set(top10["canonical_coordinate"]),
            key=lambda item: (
                int(
                    dataset_frame[
                        dataset_frame["canonical_coordinate"] == item
                    ]["rank"].min()
                ),
                canonical_sort_key(str(item)),
            ),
        )
        union_order = {item: index for index, item in enumerate(union_items)}
        subset = dataset_frame[
            dataset_frame["canonical_coordinate"].isin(union_items)
        ].copy()
        subset["union_order"] = subset["canonical_coordinate"].map(union_order)
        subset["in_top10"] = subset["rank"] <= 10
        heatmap_rows.extend(
            subset[
                [
                    "dataset",
                    "union_order",
                    "canonical_coordinate",
                    "chain",
                    "configured_group_label",
                    "run",
                    "run_index",
                    "rank",
                    "effect",
                    "sign",
                    "absolute_effect",
                    "in_top10",
                ]
            ].to_dict("records")
        )
    heatmap = pd.DataFrame(heatmap_rows).sort_values(
        ["dataset", "union_order", "run_index"], kind="stable"
    )
    return (
        top1.reset_index(drop=True),
        jaccard.reset_index(drop=True),
        heatmap.reset_index(drop=True),
        tie_diagnostics,
    )


def prepare_configured_group_supports(
    replicates: pd.DataFrame,
    *,
    replicate_source: Path,
) -> pd.DataFrame:
    """Prepare shared support-size metadata for Figures 3 and 6."""

    items = _base_item_annotations(replicates, replicate_source)
    rows: list[dict[str, object]] = []
    for (dataset, chain, label, region_id), group in items.groupby(
        ["dataset", "chain", "domain", "region_id"], sort=True
    ):
        members = tuple(
            sorted(group["canonical_item_key"].astype(str), key=canonical_sort_key)
        )
        rows.append(
            {
                "dataset": str(dataset),
                "chain": str(chain),
                "configured_group_label": str(label),
                "region_id": str(region_id),
                "support_size": len(members),
                "support_membership_hash": membership_hash(members),
                "support_membership": _represent_members(members),
                "is_unassigned_remainder": str(label).strip().lower()
                == "unassigned",
                "is_configured_catalytic_group": str(label)
                in CATALYTIC_GROUPS,
            }
        )
    result = pd.DataFrame(rows).sort_values(
        ["dataset", "chain", "configured_group_label"], kind="stable"
    )
    if len(result) != 32:
        raise ValueError(f"expected 32 configured-group rows, got {len(result)}")
    return result.reset_index(drop=True)


def _verify_v9b(paths: SourcePaths) -> dict[str, str]:
    manifest = _load_json(paths.v9b_manifest)
    recorded = manifest.get("artifact_sha256")
    if not isinstance(recorded, dict):
        raise ValueError("V9B manifest lacks artifact_sha256")
    files = {
        "v9b_partial_presence_fixture.csv": paths.v9b_fixture,
        "v9b_partial_presence_expected.csv": paths.v9b_expected,
        "v9b_partial_presence_results.json": paths.v9b_results,
    }
    hashes: dict[str, str] = {}
    for key, path in files.items():
        actual = sha256_file(path)
        hashes[key] = actual
        if recorded.get(key) != actual:
            raise ValueError(
                f"V9B hash mismatch for {key}: {recorded.get(key)} != {actual}"
            )
    results = _load_json(paths.v9b_results)
    findings = results.get("semantic_findings")
    if not isinstance(findings, dict):
        raise ValueError("V9B results lack semantic_findings")
    return hashes


def _write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(
        path,
        index=False,
        lineterminator="\n",
        float_format="%.17g",
    )


def _source_record(path: Path, role: str, base: Path | None = None) -> dict[str, object]:
    display = str(path.relative_to(base)) if base is not None else str(path)
    return {
        "role": role,
        "path": display,
        "sha256": sha256_file(path),
        "bytes": path.stat().st_size,
    }


def prepare_all(
    sources: SourcePaths,
    output_dir: Path,
    *,
    script_path: Path | None = None,
) -> dict[str, Path]:
    """Build every Phase 3E table and the deterministic manifest."""

    missing = [path for path in sources.all_inputs() if not path.is_file()]
    if missing:
        raise FileNotFoundError(f"missing authoritative inputs: {missing}")
    _verify_v9b(sources)

    stride = pd.read_parquet(sources.stride_table)
    replicates = pd.read_parquet(sources.replicate_table)
    residue = pd.read_parquet(sources.residue_landscape)
    positions = pd.read_parquet(sources.position_conservation)

    unique, paired, unique_diagnostics = prepare_unique_supports(
        stride,
        replicates,
        stride_source=sources.stride_table,
        replicate_source=sources.replicate_table,
    )
    observed_crossing = prepare_observed_crossing_sensitivity(unique)
    (
        v9a_support,
        v9a_association,
        v9a_aggregate,
        v9a_crossing,
        v9a_alias,
    ) = prepare_v9a_tables(sources)
    top1, jaccard, heatmap, tie_diagnostics = prepare_figure4_tables(
        replicates, replicate_source=sources.replicate_table
    )
    status, coordinate_summary, row_order = prepare_figure5_tables(
        residue,
        positions,
        residue_source=sources.residue_landscape,
        position_source=sources.position_conservation,
    )
    group_supports = prepare_configured_group_supports(
        replicates, replicate_source=sources.replicate_table
    )

    outputs: dict[str, pd.DataFrame] = {
        "figure2_unique_supports.csv": unique,
        "figure2_item_first_aggregate.csv": paired,
        "figure2_observed_first_crossing_sensitivity.csv": observed_crossing,
        "figure2_v9a_support_size.csv": v9a_support,
        "figure2_v9a_support_size_association.csv": v9a_association,
        "figure2_v9a_item_first_aggregate.csv": v9a_aggregate,
        "figure2_v9a_first_crossing.csv": v9a_crossing,
        "figure2_v9a_alias_provenance.csv": v9a_alias,
        "figure4_top1.csv": top1,
        "figure4_topn_jaccard.csv": jaccard,
        "figure4_top10_heatmap.csv": heatmap,
        "figure5_status_matrix.csv": status,
        "figure5_coordinate_summary.csv": coordinate_summary,
        "figure5_row_order.csv": row_order,
        "configured_group_supports.csv": group_supports,
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    written: dict[str, Path] = {}
    for filename, frame in outputs.items():
        path = output_dir / filename
        _write_csv(frame, path)
        written[filename] = path

    script = script_path or Path(__file__).resolve()
    frozen_records = [
        _source_record(sources.stride_table, "observed_support_profiles", sources.frozen_root),
        _source_record(sources.replicate_table, "per_run_effects_and_memberships", sources.frozen_root),
        _source_record(sources.residue_landscape, "observed_item_rho", sources.frozen_root),
        _source_record(sources.position_conservation, "cross_dataset_coordinate_union", sources.frozen_root),
    ]
    v9a_records = [
        _source_record(path, "completed_v9a_artifact", sources.v9a_dir)
        for path in (
            sources.v9a_snapshot,
            sources.v9a_results,
            sources.v9a_support_size,
            sources.v9a_first_crossing,
            sources.v9a_alias_effects,
            sources.v9a_manifest,
        )
    ]
    v9b_records = [
        _source_record(path, "completed_v9b_artifact", sources.v9b_dir)
        for path in (
            sources.v9b_fixture,
            sources.v9b_expected,
            sources.v9b_results,
            sources.v9b_manifest,
        )
    ]
    manifest = {
        "schema_version": "phase3e.presentation_data.1",
        "purpose": "deterministic_presentation_only_data_preparation",
        "source_commits": {
            "stride": STRIDE_COMMIT,
            "stride_downstream": DOWNSTREAM_COMMIT,
            "presentation_base": PRESENTATION_BASE_COMMIT,
        },
        "environment": {
            "python_version": platform.python_version(),
            "pandas_version": pd.__version__,
            "numpy_version": np.__version__,
        },
        "script": {
            "path": script.name,
            "sha256": sha256_file(script),
        },
        "cutoffs": list(CUTOFFS),
        "ranking": {
            "effect_column": "r",
            "primary_order": "absolute_effect_descending",
            "tie_break": "canonical_coordinate_ascending",
            **tie_diagnostics,
        },
        "figure5_row_order": {
            "n_present": "descending",
            "n_crossing_at_0.5": "descending",
            "chain": "ascending",
            "canonical_coordinate_number": "ascending",
            "canonical_coordinate": "ascending_final_tie_break",
        },
        "unique_support_diagnostics": unique_diagnostics,
        "sources": {
            "frozen": frozen_records,
            "v9a": v9a_records,
            "v9b": v9b_records,
        },
        "outputs": {
            filename: {
                "rows": int(len(outputs[filename])),
                "sha256": sha256_file(path),
                "sort_keys": list(OUTPUT_SORT_KEYS[filename]),
            }
            for filename, path in sorted(written.items())
        },
    }
    manifest_path = output_dir / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    written["manifest.json"] = manifest_path
    return written


def default_sources(root: Path) -> SourcePaths:
    configured = resolve_external_sources(outputs_root=root)
    return SourcePaths(
        frozen_root=configured.outputs_root,
        v9a_dir=configured.v9a_artifacts,
        v9b_dir=configured.v9b_artifacts,
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Prepare Phase 3E presentation-only data tables."
    )
    parser.add_argument("--root", "--outputs-root", dest="outputs_root", type=Path)
    parser.add_argument("--production-root", type=Path)
    parser.add_argument("--downstream-root", type=Path)
    parser.add_argument("--manuscript-path", type=Path)
    parser.add_argument("--v9a-dir", "--v9a-artifacts", dest="v9a_artifacts", type=Path)
    parser.add_argument("--v9b-dir", "--v9b-artifacts", dest="v9b_artifacts", type=Path)
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path(__file__).resolve().parent / "derived_data" / "phase3e",
    )
    args = parser.parse_args()
    configured = resolve_external_sources(
        production_root=args.production_root,
        downstream_root=args.downstream_root,
        outputs_root=args.outputs_root,
        manuscript_path=args.manuscript_path,
        v9a_artifacts=args.v9a_artifacts,
        v9b_artifacts=args.v9b_artifacts,
    )
    sources = SourcePaths(
        frozen_root=configured.outputs_root,
        v9a_dir=configured.v9a_artifacts,
        v9b_dir=configured.v9b_artifacts,
    )
    written = prepare_all(sources, args.output_dir)
    for filename, path in sorted(written.items()):
        print(f"{filename}: {path}")


if __name__ == "__main__":
    main()
