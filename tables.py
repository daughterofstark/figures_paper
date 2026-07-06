"""Manuscript tables T1–T5 — export the S7 tables as LaTeX / HTML / Markdown.

Reads the S7 manuscript tables verbatim (no recomputation), selects a
manuscript-friendly column subset with readable headers, rounds floats for
display, and writes booktabs LaTeX, HTML, and Markdown via
:func:`utils.export_table`.

Consumes: ``outputs_s7/T1…T5_*.parquet``.
"""
from __future__ import annotations

from pathlib import Path

import pandas as pd

from figure_config import Paths
from utils import export_table

# (s7 stem, output stem, caption, label, [columns or None for all])
_TABLES = [
    (
        "T1_per_serotype_summary", "table1_per_serotype_summary",
        "Per-serotype reproducibility summary at the provisional gate.",
        "tab:per_serotype",
        ["serotype", "n_loci", "n_gated_residue", "n_mechanisms", "n_signed",
         "n_mixed", "n_signed_significant", "frac_mixed", "rho_residue_median",
         "rho_residue_q1", "rho_residue_q3"],
    ),
    (
        "T2_domain_rho_signed_effect", "table2_domain_rho_signed_effect",
        "Domain-scale reproducibility and signed effect per serotype.",
        "tab:domain_effect",
        ["serotype", "chain", "domain", "rho_domain", "coherence_domain",
         "is_coherent", "beta_weighted_mean", "beta_weighted_se", "n_signed",
         "n_significant_fdr"],
    ),
    (
        "T3_catalytic_cross_serotype", "table3_catalytic_cross_serotype",
        "Catalytic-machinery reproducibility across serotypes.",
        "tab:catalytic",
        ["domain", "serotype", "chain", "rho_domain", "beta_domain"],
    ),
    (
        "T4_top_shared_signed_positions", "table4_top_shared_signed_positions",
        "Top reproducible signed positions shared across serotypes.",
        "tab:shared_positions",
        ["canon_label", "chain", "domain", "n_serotypes_reproducible",
         "n_serotypes_signed_reproducible", "frac_reproducible",
         "conservation_class", "is_catalytic_triad", "rho_residue_median"],
    ),
    (
        "T5_variance_component_budget", "table5_variance_component_budget",
        "Per-domain variance-component budget (τ² vs σ̄²).",
        "tab:variance_budget",
        ["serotype", "chain", "domain", "tau2", "sigma2_bar", "frac_tau2",
         "frac_sigma2", "tau2_sigma2_ratio", "variance_regime"],
    ),
]


def _round_floats(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    for c in out.columns:
        if pd.api.types.is_float_dtype(out[c]):
            out[c] = out[c].round(3)
    return out


def build(paths: Paths) -> list[Path]:
    written: list[Path] = []
    for s7_stem, out_stem, caption, label, cols in _TABLES:
        src = paths.s7(f"{s7_stem}.parquet")
        if not src.is_file():
            raise FileNotFoundError(f"missing S7 table: {src}")
        df = pd.read_parquet(src)
        if cols:
            df = df[[c for c in cols if c in df.columns]]
        df = _round_floats(df)
        written += export_table(df, paths, out_stem, caption=caption, label=label)
    return written


if __name__ == "__main__":
    from figure_config import resolve_paths

    for p in build(resolve_paths()):
        print(p)
