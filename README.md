# paper_figures

Publication-quality figure suite for the STRIDE dengue NS2B–NS3 reproducibility
analysis. This project is **presentation-only**: it reads the frozen pipeline's
parquet outputs (`outputs/`, `outputs_s1a/ … outputs_s7/`) and renders figures and
tables. It never recomputes statistics, alters thresholds, or touches raw STRIDE data.

## The manuscript at a glance (6 main + 2 supplementary)

The suite is designed as one coherent visual narrative: *where* reproducible dynamics
live → *at what scale* they resolve → *which module* carries them → *what the signed
effects are* → *how they behave across serotypes* → *what limits them*.

| Figure | One-sentence biological message | Built from |
| --- | --- | --- |
| **Fig 1** — reproducibility landscape | Residue-level reproducibility is heterogeneous and distributed across NS2B–NS3. | `s5/position_conservation` |
| **Fig 2** — spatial resolution | Catalytic-region ρ plateaus after residue aggregation, while gates remain mostly residue/SS. | `s7/F7` (residue→SS→motif→domain) + `s7/F2` |
| **Fig 3** — domain/coherence audit | Catalytic domains are directionally coherent but are not the highest-ρ regions. | `s7/F3` + `s7/F8` |
| **Fig 4** — signed-effect volcano | FDR-significant signed effects are widespread across NS2B–NS3, with catalytic examples highlighted. | `s4/significance_screen` |
| **Fig 5** — cross-serotype audit | Cross-serotype reproducibility is heterogeneous; positions conserved in more serotypes have higher median ρ. | `s5/position_conservation` + `s7/F1` |
| **Fig 6** — variance decomposition | Reproducibility is limited by between-replicate variance (τ²), not sampling. | `s7/F6` |
| **Supp S1** — per-serotype landscape | Per-serotype residue-level ρ landscapes are heterogeneous, with catalytic residues as landmarks. | `s7/F1` |
| **Supp S2** — signed-β Manhattan | Signed effects are distributed along the sequence, complementing the Fig 4 volcano. | `s4/significance_screen` |

Each module's docstring records its biological message, rationale, source tables, and
Main-vs-Supplementary justification for future maintenance.

## Global visual system

**One colour = one meaning** (see `styles.py`): catalytic = vermillion (reserved);
serotype identity = blue/sky/green/purple; reproducibility ρ = viridis; conservation =
blue ramp; correlation = cividis; FDR significance = darkness (ink vs light grey);
variance τ² = orange, σ̄² = neutral grey. Colours are never reused for unrelated
concepts.

**One label rule** (`utils.objective_labels`): for signed-effect plots, label an
observation only when it is FDR-significant and its |score| is in the top decile, capped
at 4 by |score| and de-overlapped with adjustText. Used in Figure 4 and Supplementary
Figure S2.

**Typography / lettering**: a single rcParams style (`styles.apply_style`) fixes font
family, sizes, tick sizing and panel-letter placement (`utils.panel_letter`) across
every figure. Declarative, data-justified titles throughout.

## Layout

```
paper_figures/
  README.md          CHANGELOG.md       figure_config.py   styles.py   utils.py
  build_all.py       tables.py          conftest.py
  figure1.py … figure6.py               # 6 main figures
  supplementary_figure1.py, supplementary_figure2.py
  requirements.txt   requirements-dev.txt
  tests/             # regression tests + real-like multi-chain fixture
  output/            # generated artifacts (git-ignored)
```

## Build

```
python build_all.py                       # uses ./outputs* or $STRIDE_OUTPUTS_ROOT
python build_all.py --root /path/to/outs --output-dir output
python figure5.py                         # any module builds standalone
```

For **every figure**: `output/figureN_<name>.{svg,pdf,png}` (main) and
`output/Supplementary_Figure_S{1,2}.{svg,pdf,png}` (supplementary) — vector SVG + PDF
for typesetting and 600 dpi PNG for viewing. Tables T1–T5 export as `.tex/.html/.md`.

## Reproducibility

Deterministic and byte-reproducible: metadata/timestamps are stripped on save,
adjustText runs with a fixed iteration cap, and clustering (Fig 5B) is deterministic.
Two builds of the full suite are identical (0/39 files differ). The suite builds
cleanly under `PYTHONWARNINGS=error` on both the multi-chain real-like outputs and the
single-chain synthetic example (figures degrade gracefully when a chain/domain is
absent).
