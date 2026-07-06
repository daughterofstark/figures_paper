# paper_figures

Publication-quality figures and manuscript tables for the STRIDE dengue
downstream analysis — built **entirely from the pipeline's existing outputs**.

This is a **presentation-only** project. It is deliberately separate from the
`stride-dengue-analysis` framework (S0–S7), which remains the authoritative,
reproducible computational pipeline. `paper_figures` consumes the tables that
pipeline already wrote and redesigns them for manuscript submission.

> It never reads raw STRIDE outputs, never reruns S0–S7, never recomputes a
> statistic, and never performs new inference. Every number shown is read verbatim
> from a pipeline output table.

---

## Required inputs

A completed STRIDE downstream run — i.e. a directory (the *outputs root*) that
contains:

```
<root>/
  outputs_s2/   outputs_s3/   outputs_s4/
  outputs_s5/   outputs_s6/   outputs_s7/
```

The project reads primarily from `outputs_s7/` (the S7 reporting layer's
prepared per-figure tables `F1…F8_*.parquet` and manuscript tables
`T1…T5_*.parquet`), and enriches two figures from upstream tables that carry extra
columns:

| Figure | Also reads |
| --- | --- |
| Figure 4 (volcano) | `outputs_s4/significance_screen.parquet` (β + BH-adjusted p) |
| Figure 5 (conservation) | `outputs_s5/position_conservation.parquet` (per-class + catalytic flags) |

Run the pipeline through S7 first (see the framework's README). Then point this
project at the outputs root (see **Regenerate** below).

## Installation

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
```

Dependencies: `matplotlib`, `pandas`, `numpy`, `pyarrow` (parquet), `scipy`, and
`adjustText` (non-overlapping labels; the figures degrade gracefully without it).
`seaborn` is intentionally not used.

## Regenerate all figures and tables

```bash
python build_all.py --root /path/to/outputs-root
# or, if the analysis repo sits next to this project, just:
python build_all.py
# or set the environment variable:
STRIDE_OUTPUTS_ROOT=/path/to/outputs-root python build_all.py
```

Every artifact is written to `output/`. Individual figures can be built alone,
e.g. `python figure4.py`.

## Generated outputs

For **every figure**: `output/figureN_<name>.{svg,pdf,png}` — vector SVG + PDF for
typesetting and a 600 dpi PNG for quick viewing. The five **supplementary**
figures follow the same format as `output/Supplementary_Figure_S{1..5}.{svg,pdf,png}`.

For **every table**: `output/tableN_<name>.{tex,html,md}` — booktabs LaTeX, HTML,
and Markdown.

## Figure ↔ manuscript-panel correspondence

Each figure redesigns the identically-numbered S7 figure (F1–F8); each table
corresponds to the S7 manuscript table (T1–T5).

| File | Redesigns | Shows | Redesign |
| --- | --- | --- | --- |
| `figure1_reproducibility_landscape` | F1 | per-residue ρ along the sequence | serotype-faceted small multiples + rolling mean + catalytic highlight (was a spaghetti plot) |
| `figure2_resolution_census` | F2 | gated loci per spatial scale, per serotype | improved stacked bars (scale-keyed sequential colour, in-bar counts) |
| `figure3_domain_serotype_heatmap` | F3 | ρ over domain × serotype | heatmap ordered by mean ρ, catalytic rows flagged, annotated cells |
| `figure4_signed_effect_volcano` | F4 | signed effect vs significance | **volcano** (β vs −log10 BH-p), FDR-accented, catalytic labelled (was an unreadable forest of CIs) |
| `figure5_cross_serotype_conservation` | F5 | conservation of reproducible positions | ranked dot plot + per-class summary bar (was unreadable) |
| `figure6_variance_composition` | F6 | τ² vs σ̄² per domain | serotype-faceted, domains grouped catalytic-first |
| `figure7_rho_vs_scale_catalytic` | F7 | ρ across spatial scale, catalytic regions | serotype-faceted trajectories + bold mean, de-cluttered |
| `figure8_coherence_vs_rho` | F8 | coherence vs ρ at domain scale | labelled catalytic domains, provisional quadrant guides |
| `table1_per_serotype_summary` | T1 | per-serotype summary | booktabs / HTML / Markdown |
| `table2_domain_rho_signed_effect` | T2 | domain ρ + signed effect | " |
| `table3_catalytic_cross_serotype` | T3 | catalytic cross-serotype | " |
| `table4_top_shared_signed_positions` | T4 | top shared signed positions | " |
| `table5_variance_component_budget` | T5 | variance-component budget | " |

### Supplementary figures (S1–S5)

New supplementary figures that only *visualize* existing outputs (Figures 1–8 are
unchanged). Filenames: `output/Supplementary_Figure_S{1..5}.{svg,pdf,png}`.

| File | Shows | Built from (existing outputs) |
| --- | --- | --- |
| `supplementary_figure1` (S1) | structural domain map of NS2B–NS3 with per-position ρ overlaid as lollipops; chains shaded, domains blocked, catalytic ringed | `outputs_s5/position_conservation` (chain, domain, `rho_residue_median`, `is_catalytic_triad`; position parsed from `canon_label`) |
| `supplementary_figure2` (S2) | conservation vs signed effect scatter (colour = domain, size = ρ), catalytic ringed, deterministic LOWESS trend | `outputs_s5/position_conservation` + `outputs_s4/significance_screen`, joined on `canon_label` |
| `supplementary_figure3` (S3) | ranked domain summary: median/mean ρ bars + coherence and τ²-fraction dots; catalytic flagged | `outputs_s7/T2_domain_rho_signed_effect` + `outputs_s7/T5_variance_component_budget` (aggregated per (chain, domain)) |
| `supplementary_figure4` (S4) | Manhattan/lollipop of signed β along the sequence; sign-coloured, size = FDR, top-10 ± labelled | `outputs_s4/significance_screen` |
| `supplementary_figure5` (S5) | 4×4 cross-serotype similarity of ρ profiles (Pearson r) with hierarchical clustering + dendrogram | `outputs_s7/F1_reproducibility_landscape` (pivoted to positions × serotypes) |

Supplementary figures are strictly presentation: S1–S4 use only pre-existing
quantities and deterministic transformations (min/max spans, group aggregation);
the LOWESS curve (S2) and the pairwise correlation + clustering (S5) are
visualization aids that introduce no new claimed statistics, per the brief. All
degrade gracefully on the single-chain synthetic example and render fully on the
real multi-chain (NS2B–NS3) outputs.

## Design conventions

- **Colour-blind-safe** Okabe–Ito qualitative palette; perceptually-uniform
  `viridis`/`cividis` for continuous fields.
- **Consistent typography** and a single type scale across every panel;
  no top/right spines, no grid clutter.
- **Panel lettering** (A, B, …) on multi-panel figures.
- **Vector-first** output (SVG/PDF); text kept as text (`svg.fonttype = none`,
  `pdf.fonttype = 42`) so labels stay editable and selectable.
- **Deterministic**: fixed SVG hash-salt and stripped timestamps mean rebuilding
  yields byte-identical SVG, PDF, and PNG.
- **Uncalibrated honesty**: any ρ*/coherence guide lines are drawn only as
  provisional references and annotated as such — the figures make no calibrated
  pass/fail claim, mirroring the pipeline.

## What this project must not do

Never edit the S0–S7 framework, its CLIs, tests, or outputs. Never recompute a
statistic or run a new hypothesis test. Improve presentation only.

## Layout

```
paper_figures/
  README.md          figure_config.py   styles.py   utils.py
  build_all.py       tables.py          conftest.py
  figure1.py … figure8.py
  supplementary_figure1.py … supplementary_figure5.py
  requirements.txt   requirements-dev.txt
  tests/             # regression tests (incl. real-data multi-chain fixture)
  output/            # generated artifacts (git-ignored)
```
