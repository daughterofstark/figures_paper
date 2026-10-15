# STRIDE paper figures

This repository contains the presentation layer for the STRIDE descriptive
methods paper. It reads frozen downstream outputs and committed Phase 3E tables,
then renders figures, tables, captions, and numeric audit records. It does not
fit the STRIDE estimator, change thresholds, or recompute scientific results.

## Paper architecture

| Display item | Descriptive question |
| --- | --- |
| Figure 1 | How does item-level magnitude stability vary across coordinates and datasets, and where are coordinates absent? |
| Figure 2 | How does support-specific magnitude stability change between an item and its first unique aggregate, and how does V9A characterize support-size behavior? |
| Figure 3 | How do configured-group magnitude rho and profile coherence differ? |
| Figure 4 | How much do single-run top-item and top-N rankings vary among run labels? |
| Figure 5 | How do coordinate presence and provisional item-level crossing vary across datasets and reference values? |
| Figure 6 | How does the fitted K=3 model divide variance between the fitted between-run component and propagated within-run uncertainty? |
| Supplementary Figure S1 | What are the complete per-dataset item-rho profiles? |
| Supplementary Figure S2 | What conditional selected-support signed associations are present under the implemented provisional screen? |

The figures use continuous rho as a support-specific magnitude-stability
summary. The value `0.5` is a provisional operational reference. Larger
supports define different L2 aggregate objects, so cross-support rho changes do
not form an evidence ladder. Signed and BH outputs are conditional descriptive
summaries. Cross-dataset recurrence distinguishes absence from presence below a
reference value.

## Inputs

External source discovery is centralized in `figure_config.py`. A command-line
argument overrides an environment variable, which overrides a documented
relative default.

Supported environment variables:

- `STRIDE_PRODUCTION_ROOT`
- `STRIDE_DOWNSTREAM_ROOT`
- `STRIDE_OUTPUTS_ROOT`
- `STRIDE_MANUSCRIPT_PATH`
- `STRIDE_V9A_ARTIFACTS`
- `STRIDE_V9B_ARTIFACTS`

Inspect the resolved paths:

```bash
python figure_config.py
```

The prepared frozen-input archive is not stored in this repository. Integration
tests skip with an explicit reason when their external inputs are unavailable.
A skipped integration test is not evidence that the external-data workflow
passed.

## Installation

```bash
git clone https://github.com/daughterofstark/figures_paper.git
cd figures_paper
python -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[dev]"
```

## Prepare Phase 3E tables

Phase 3E is the read-only presentation-data layer used by Figures 2, 4, and 5.
Run it with explicit roots:

```bash
python prepare_phase3e_data.py \
  --production-root /path/to/STRIDE \
  --downstream-root /path/to/STRIDE_downstream \
  --outputs-root /path/to/strideouts \
  --manuscript-path /path/to/STRIDE_MANUSCRIPT_DESCRIPTIVE_METHODS_WORKING.md
```

The command writes presentation tables and their manifest. It does not overwrite
the frozen STRIDE or downstream outputs.

## Build

Build the full presentation set:

```bash
python build_all.py --root /path/to/strideouts --output-dir output
```

Individual modules also build independently:

```bash
python figure3.py
python supplementary_figure2.py
python tables.py
```

Figures are exported as PDF, SVG, and high-resolution PNG. Table outputs include
machine-readable and manuscript-ready forms where supported. Plotting modules
validate required columns and denominators before rendering.

## Tests

Self-contained checks require no archived scientific data:

```bash
python -m compileall -q .
ruff check --select E9,F63,F7,F82 .
MPLBACKEND=Agg pytest -m "not integration" -q
```

Run the complete suite by supplying the external roots:

```bash
export STRIDE_PRODUCTION_ROOT=/path/to/STRIDE
export STRIDE_DOWNSTREAM_ROOT=/path/to/STRIDE_downstream
export STRIDE_OUTPUTS_ROOT=/path/to/strideouts
export STRIDE_MANUSCRIPT_PATH=/path/to/STRIDE_MANUSCRIPT_DESCRIPTIVE_METHODS_WORKING.md
export STRIDE_V9A_ARTIFACTS=/path/to/STRIDE/validation/artifacts/v9a
export STRIDE_V9B_ARTIFACTS=/path/to/STRIDE_downstream/validation_artifacts
MPLBACKEND=Agg pytest -q
```

GitHub Actions runs the self-contained tests on Python 3.10, 3.11, and 3.12. It
also verifies that absent integration data produce explicit skips.

## Reproducibility boundary

This repository verifies presentation transformations, deterministic ordering,
and agreement with frozen source values. It does not establish threshold
calibration, post-selection inference, biological homology, or mechanistic
interpretation. Those limits are part of the manuscript and should remain
visible in any reuse of the figures or tables.

## Historical audit material

Directories named `manuscript_number_audit` and `manuscript_placeholder_audit`
preserve earlier review records. They refer to a retired manuscript and figure
architecture and are not sources for the current scientific wording. Use the
active scripts, `output/` audit records, Phase 3 reports, and the current
descriptive manuscript when checking present claims.

## License

MIT. See `LICENSE`.
