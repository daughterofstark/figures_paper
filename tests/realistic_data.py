"""Generate real-like STRIDE downstream outputs for regression tests.

Unlike the tiny committed example (whose domains all sit on a single chain), the
real NS2B–NS3 protease is a **two-chain complex**, so the same domain label can
appear on more than one chain within a serotype. This generator reproduces that
condition: the domain ``"Interface"`` occurs on both ``NS2B`` and ``NS3``. That is
exactly what makes ``domain`` non-unique within a serotype and previously crashed
Figure 6 (``reindex`` on duplicate labels) and silently averaged Figure 3.

All values are deterministic functions of indices — no randomness, no real data.
Column sets mirror the true S4/S5/S7 schemas the figures consume.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

SEROTYPES = ("DENV1", "DENV2", "DENV3", "DENV4")

# (chain, domain, is_catalytic) — note "Interface" recurs on NS2B *and* NS3.
REGIONS = (
    ("NS3", "Catalytic Triad", True),
    ("NS3", "Oxyanion Loop", True),
    ("NS3", "Helicase Lobe 1", False),
    ("NS3", "Interface", False),
    ("NS2B", "Cofactor Core", False),
    ("NS2B", "Interface", False),
)

# (canon_label, chain, domain)
RESIDUES = (
    ("NS3:51", "NS3", "Catalytic Triad"),
    ("NS3:75", "NS3", "Catalytic Triad"),
    ("NS3:135", "NS3", "Catalytic Triad"),
    ("NS3:38", "NS3", "Oxyanion Loop"),
    ("NS3:200", "NS3", "Helicase Lobe 1"),
    ("NS3:150", "NS3", "Interface"),
    ("NS2B:80", "NS2B", "Cofactor Core"),
    ("NS2B:52", "NS2B", "Interface"),
)

SCALES = (
    (0, "residue"),
    (1, "secondary_structure"),
    (2, "motif"),
    (3, "domain"),
    (4, "chain"),
    (5, "protein"),
    (6, "complex"),
)
_CATALYTIC = {"Catalytic Triad", "Oxyanion Loop"}
_TIER = "A_licensed"


def _sero_idx(s: str) -> int:
    return SEROTYPES.index(s)


def _rho_domain(s: str, chain: str, dom: str) -> float:
    base = 0.9 if dom in _CATALYTIC else 0.55
    return round(base - 0.03 * _sero_idx(s) - (0.05 if chain == "NS2B" else 0.0), 3)


# --------------------------------------------------------------------------
# Per-figure frames (S7 prepared data)
# --------------------------------------------------------------------------
def _f1() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        si = _sero_idx(s)
        for ri, (canon, chain, dom) in enumerate(RESIDUES):
            # base reproducibility plus a serotype×residue interaction so the
            # per-serotype ρ profiles have genuinely different *shapes* (not just
            # a constant offset) — makes the S5 correlation matrix informative.
            base = _rho_domain(s, chain, dom) + 0.01
            wobble = 0.06 * np.sin(0.9 * ri + 1.3 * si) + 0.03 * ((ri * (si + 1)) % 3 - 1)
            rho = min(0.99, max(0.05, base + wobble))
            rows.append(dict(serotype=s, canon_label=canon, chain=chain,
                             domain=dom, rho_residue=round(rho, 3), tier=_TIER))
    return pd.DataFrame(rows)


def _f2() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for idx, (si, sl) in enumerate([(0, "residue"), (3, "domain")]):
            rows.append(dict(serotype=s, gated_scale_level=sl,
                             gated_scale_index=si, n_loci=3 + idx + _sero_idx(s)))
    return pd.DataFrame(rows)


def _f3() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for chain, dom, cat in REGIONS:
            rows.append(dict(domain=dom, serotype=s, chain=chain,
                             rho_domain=_rho_domain(s, chain, dom),
                             is_catalytic_domain=cat, tier=_TIER))
    return pd.DataFrame(rows)


def _f4() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for canon, chain, dom in RESIDUES:
            if dom in _CATALYTIC or chain == "NS2B":
                b = 0.3 - 0.05 * _sero_idx(s)
                rows.append(dict(serotype=s, canon_label=canon, domain=dom,
                                 beta_signed=round(b, 3),
                                 beta_ci_lower=round(b - 0.1, 3),
                                 beta_ci_upper=round(b + 0.1, 3),
                                 significant_fdr=dom in _CATALYTIC, tier=_TIER))
    return pd.DataFrame(rows)


def _f5() -> pd.DataFrame:
    rows = []
    for canon, chain, dom in RESIDUES:
        n = 4 if dom in _CATALYTIC else 2
        rows.append(dict(canon_label=canon, chain=chain, domain=dom,
                         n_serotypes_reproducible=n,
                         frac_reproducible=round(n / 4, 3),
                         conservation_class=("reproducible_all" if n == 4
                                             else "reproducible_some"),
                         is_serotype_divergent=(n < 4), tier=_TIER))
    return pd.DataFrame(rows)


def _f6() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for chain, dom, cat in REGIONS:
            ft = 0.4 + 0.05 * _sero_idx(s) + (0.1 if chain == "NS2B" else 0.0)
            ft = min(ft, 0.9)
            rows.append(dict(serotype=s, chain=chain, domain=dom,
                             frac_tau2=round(ft, 3), frac_sigma2=round(1 - ft, 3),
                             variance_regime=("replicate_limited" if ft > 0.5
                                              else "sampling_limited"),
                             tier=_TIER))
    return pd.DataFrame(rows)


def _f7() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for canon, chain, dom in RESIDUES:
            if dom not in _CATALYTIC:
                continue
            for si, sl in SCALES:
                rho = 0.7 + 0.03 * si - 0.02 * _sero_idx(s)
                rows.append(dict(serotype=s, canon_label=canon, domain=dom,
                                 scale_index=si, scale_level=sl,
                                 rho=round(min(rho, 0.99), 3), tier=_TIER))
    return pd.DataFrame(rows)


def _f8() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for chain, dom, cat in REGIONS:
            rho = _rho_domain(s, chain, dom)
            coh = round(min(0.95, rho + 0.05), 3)
            rows.append(dict(serotype=s, chain=chain, domain=dom,
                             rho_domain=rho, coherence_domain=coh,
                             is_coherent=coh >= 0.6, tier=_TIER))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------
# Manuscript tables (T1–T5) and the two upstream enrichment tables
# --------------------------------------------------------------------------
def _t1() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        rows.append(dict(serotype=s, n_loci=8, n_gated_residue=6, n_mechanisms=5,
                         n_signed=4, n_mixed=1, n_signed_significant=3,
                         frac_mixed=0.2, rho_residue_median=0.8, rho_residue_q1=0.6,
                         rho_residue_q3=0.9, rho_residue_min=0.4,
                         rho_residue_max=0.95, rho_star=0.5))
    return pd.DataFrame(rows)


def _t2() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for chain, dom, cat in REGIONS:
            rho = _rho_domain(s, chain, dom)
            rows.append(dict(serotype=s, chain=chain, domain=dom, rho_domain=rho,
                             beta_domain=0.3, coherence_domain=min(0.95, rho + 0.05),
                             is_coherent=rho + 0.05 >= 0.6, beta_weighted_mean=0.28,
                             beta_weighted_se=0.05, n_signed=2, n_significant_fdr=1,
                             tier=_TIER))
    return pd.DataFrame(rows)


def _t3() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for chain, dom, cat in REGIONS:
            if not cat:
                continue
            rows.append(dict(domain=dom, serotype=s, chain=chain,
                             rho_domain=_rho_domain(s, chain, dom),
                             beta_domain=0.3, is_catalytic_domain=True, tier=_TIER))
    return pd.DataFrame(rows)


def _t4() -> pd.DataFrame:
    rows = []
    for canon, chain, dom in RESIDUES:
        cat = dom in _CATALYTIC
        n = 4 if cat else 2
        rows.append(dict(canon_label=canon, chain=chain, domain=dom,
                         n_serotypes_present=4, n_serotypes_reproducible=n,
                         n_serotypes_signed_reproducible=(n if cat else 1),
                         frac_reproducible=round(n / 4, 3),
                         conservation_class=("reproducible_all" if n == 4
                                             else "reproducible_some"),
                         is_serotype_divergent=(n < 4), is_catalytic_triad=cat,
                         rho_residue_median=0.8, tier=_TIER))
    return pd.DataFrame(rows)


def _t5() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        for chain, dom, cat in REGIONS:
            ft = min(0.9, 0.4 + 0.05 * _sero_idx(s) + (0.1 if chain == "NS2B" else 0))
            rows.append(dict(serotype=s, chain=chain, domain=dom, tau2=0.2,
                             sigma2_bar=0.3, frac_tau2=round(ft, 3),
                             frac_sigma2=round(1 - ft, 3),
                             tau2_sigma2_ratio=round(ft / (1 - ft), 3),
                             variance_regime=("replicate_limited" if ft > 0.5
                                              else "sampling_limited"), tier=_TIER))
    return pd.DataFrame(rows)


def _significance_screen() -> pd.DataFrame:
    rows = []
    for s in SEROTYPES:
        si = _sero_idx(s)
        for ri, (canon, chain, dom) in enumerate(RESIDUES):
            signed = dom in _CATALYTIC or chain == "NS2B"
            # per-residue variation so effects are not coincident within a serotype
            b = 0.30 - 0.05 * si + 0.05 * np.sin(1.1 * ri + 0.7 * si)
            pbh = (0.003 + 0.05 * abs(np.sin(0.6 * ri + 0.9 * si))
                   if dom in _CATALYTIC else 0.30)
            rows.append(dict(serotype=s, canon_label=canon, chain=chain, domain=dom,
                             gated_scale_level="residue", tier=_TIER,
                             direction=("increase" if signed else "mixed"),
                             is_signed=signed,
                             beta_signed=(round(b, 3) if signed else None),
                             beta_se=0.05,
                             beta_ci_lower=(round(b - 0.1, 3) if signed else None),
                             beta_ci_upper=(round(b + 0.1, 3) if signed else None),
                             ci_excludes_zero=signed,
                             z_score=(round(b / 0.05, 3) if signed else 0.0),
                             p_value=(round(pbh * 0.8, 4) if dom in _CATALYTIC else 0.2),
                             p_value_bh=round(pbh, 4),
                             significant_raw=dom in _CATALYTIC,
                             significant_fdr=dom in _CATALYTIC))
    return pd.DataFrame(rows)


def _position_conservation() -> pd.DataFrame:
    rows = []
    for ri, (canon, chain, dom) in enumerate(RESIDUES):
        cat = dom in _CATALYTIC
        # catalytic residues are highly (but not always perfectly) conserved;
        # vary between 3 and 4 so they spread rather than stacking at one x.
        n = (4 if ri % 2 == 0 else 3) if cat else 2
        rows.append(dict(canon_label=canon, chain=chain, domain=dom,
                         n_serotypes_total=4, n_serotypes_present=4,
                         serotypes_present="DENV1;DENV2;DENV3;DENV4",
                         in_all_serotypes=True, n_serotypes_reproducible=n,
                         n_serotypes_signed_reproducible=(n if cat else 1),
                         frac_reproducible=round(n / 4, 3),
                         conservation_class=("reproducible_all" if n == 4
                                             else "reproducible_majority" if n == 3
                                             else "reproducible_some"),
                         is_serotype_divergent=(n < 4), is_catalytic_triad=cat,
                         rho_residue_min=0.6, rho_residue_median=0.8,
                         rho_residue_max=0.95, rho_star=0.5,
                         is_provisional_rho_star=True, tier=_TIER))
    return pd.DataFrame(rows)


_S7_FRAMES = {
    "F1_reproducibility_landscape": _f1,
    "F2_resolution_census": _f2,
    "F3_domain_serotype_rho_heatmap": _f3,
    "F4_signed_effect_forest": _f4,
    "F5_cross_serotype_conservation": _f5,
    "F6_variance_composition": _f6,
    "F7_rho_vs_scale_catalytic": _f7,
    "F8_coherence_vs_rho": _f8,
    "T1_per_serotype_summary": _t1,
    "T2_domain_rho_signed_effect": _t2,
    "T3_catalytic_cross_serotype": _t3,
    "T4_top_shared_signed_positions": _t4,
    "T5_variance_component_budget": _t5,
}


def write_realistic_outputs(root: Path) -> Path:
    """Materialise a real-like outputs root (two-chain, duplicate-domain) and return it."""
    root = Path(root)
    s7 = root / "outputs_s7"
    s7.mkdir(parents=True, exist_ok=True)
    for stem, fn in _S7_FRAMES.items():
        fn().to_parquet(s7 / f"{stem}.parquet", index=False)
    s4 = root / "outputs_s4"
    s4.mkdir(parents=True, exist_ok=True)
    _significance_screen().to_parquet(s4 / "significance_screen.parquet", index=False)
    s5 = root / "outputs_s5"
    s5.mkdir(parents=True, exist_ok=True)
    _position_conservation().to_parquet(s5 / "position_conservation.parquet", index=False)
    return root
