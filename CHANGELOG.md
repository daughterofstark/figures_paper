# Changelog — consolidation redesign

Transforms the 13-figure suite into a coherent **6 main + 2 supplementary** manuscript
by implementing the PI-level review. Purely presentational: no analysis, statistic, or
threshold was changed; only existing pipeline outputs are visualized.

## Old → New figure mapping (for updating manuscript text)

| Old | New | Action |
| --- | --- | --- |
| old Supp S1 (structural map) | **Fig 1** | Promoted to main; declarative title; median-ρ overlay. |
| old Fig 1 (per-serotype landscape) | **Supp S1** | Demoted to supplement; kept as distinct per-serotype detail. |
| old Fig 7 (ρ vs scale) | **Fig 2A** | Truncated to residue→SS→motif→domain; trivial coarse scales removed. |
| old Fig 2 (resolution census) | **Fig 2B** | Merged into the spatial-resolution figure. |
| old Fig 3 (domain ρ heatmap) | **Fig 3A** | Merged into the catalytic-core figure. |
| old Fig 8 (coherence vs ρ) | **Fig 3B** | Merged into the catalytic-core figure. |
| old Fig 4 (volcano) | **Fig 4** | New colour semantics + single global label rule; FDR line labelled, no implied gate. |
| old Fig 5 (conservation dots) | **Fig 5A** | Merged into the headline cross-serotype figure; bubble-size encoding dropped. |
| old Supp S5 (serotype matrix) | **Fig 5B** | Merged into the headline figure (clustered heatmap + dendrogram). |
| old Supp S2 (conservation vs effect) | **Fig 5C** | Merged in and promoted — the manuscript's headline result. |
| old Fig 6 (variance composition) | **Fig 6** | Kept as full τ²/σ̄² decomposition; τ² recoloured orange; declarative title. |
| old Supp S4 (signed-β Manhattan) | **Supp S2** | Kept as positional supplement; recoloured to match Fig 4 significance semantics. |
| old Supp S3 (domain-summary ranking) | — | **Deleted** (Panel A duplicated Fig 3; Panel B duplicated Figs 6 & 8). |

Deleted modules: `figure7.py`, `figure8.py`, `supplementary_figure3.py`,
`supplementary_figure4.py`, `supplementary_figure5.py`.

## Global changes

- **Colour semantics rewritten (`styles.py`)** to one-colour-one-meaning: vermillion
  reserved for catalytic only (was also significance and τ²); DENV2 moved off orange to
  sky; τ² → orange, σ̄² → neutral grey; FDR significance → darkness (ink vs light grey);
  removed the coherence colour (coherence is now an axis, not a hue).
- **Single objective label rule** added (`utils.objective_labels`): catalytic OR
  (FDR-significant AND top-decile |score|), capped at 8; replaces all per-figure manual
  label choices. Applied in Figs 4, 5C, Supp S2.
- **Shared helpers** (`utils.py`): `lowess` (deterministic visual trend), `chain_bands`
  (sequence-axis chain shading), `objective_labels`; no duplicated plotting code across
  the merged multi-panel figures.
- **Titles** rewritten as declarative, data-justified biological statements; **axes**
  decluttered (fewer ticks, fewer decimals, self-explanatory labels); **legends** merged
  to one per figure; **panel lettering** standardised via `utils.panel_letter`.
- **Build pipeline (`build_all.py`)**: `FIGURE_MODULES = figure1..6`,
  `SUPPLEMENTARY_MODULES = supplementary_figure1..2`, tables unchanged. Full suite =
  8 figures × 3 formats + 5 tables × 3 formats = 39 files.
- **Tests** updated to the new set (6 passing); the real-like fixture gained per-residue
  variation in ρ, β and p so the merged panels are exercised meaningfully.

## Invariants preserved

Deterministic, byte-reproducible output (0/39 files differ across builds); vector
SVG/PDF + 600 dpi PNG; colour-blind-safe palette; no warnings under
`PYTHONWARNINGS=error`; backward-compatible with the single-chain synthetic example;
`(chain, domain)` keying retained so duplicate domain labels across chains never crash
or get averaged together.
