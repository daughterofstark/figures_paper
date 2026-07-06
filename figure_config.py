"""Configuration for paper_figures: input locations, sizes, and constants.

The project consumes the *already-generated* STRIDE downstream outputs. It never
reads raw STRIDE files and never reruns or recomputes any stage. All input paths
are resolved from a single ``outputs root`` — the directory that contains the
``outputs_s2/`` … ``outputs_s7/`` folders — so the project can be pointed at any
completed run.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path

# --------------------------------------------------------------------------
# Figure geometry (journal column widths, in millimetres → inches)
# --------------------------------------------------------------------------
MM = 1.0 / 25.4
SINGLE_COL_MM = 89.0  # Nature single column
ONEHALF_COL_MM = 120.0
DOUBLE_COL_MM = 183.0  # Nature double column

SINGLE_COL = SINGLE_COL_MM * MM
ONEHALF_COL = ONEHALF_COL_MM * MM
DOUBLE_COL = DOUBLE_COL_MM * MM

DPI = 600
FORMATS = ("svg", "pdf", "png")

# --------------------------------------------------------------------------
# Domain / scale biology (labels only — no analysis)
# --------------------------------------------------------------------------
CATALYTIC_DOMAINS = ("Catalytic Triad", "Oxyanion Loop")
SEROTYPE_ORDER = ("DENV1", "DENV2", "DENV3", "DENV4")
SCALE_ORDER = (
    "residue",
    "secondary_structure",
    "motif",
    "domain",
    "chain",
    "protein",
    "complex",
)
SCALE_SHORT = {
    "residue": "res",
    "secondary_structure": "SS",
    "motif": "motif",
    "domain": "dom",
    "chain": "chain",
    "protein": "prot",
    "complex": "cplx",
}
CONSERVATION_ORDER = (
    "reproducible_all",
    "reproducible_majority",
    "reproducible_some",
    "reproducible_none",
)

#: Provisional (uncalibrated) reference lines carried through from the pipeline.
#: These are drawn only as visual guides and annotated as provisional; the
#: project makes no calibrated claim.
PROVISIONAL_RHO_STAR = 0.5
PROVISIONAL_COHERENCE_THRESHOLD = 0.6


@dataclass(frozen=True)
class Paths:
    """Resolved input directories and the figure output directory."""

    root: Path
    output: Path

    def stage(self, n: str) -> Path:
        """Directory for a stage's outputs, e.g. ``stage("s7")``."""
        # S0 writes to ``outputs/``; S2–S7 write to ``outputs_sN/``.
        name = "outputs" if n in ("s0", "outputs") else f"outputs_{n}"
        return self.root / name

    def s7(self, stem: str) -> Path:
        return self.stage("s7") / stem

    def table(self, stage: str, stem: str) -> Path:
        return self.stage(stage) / stem


def default_root() -> Path:
    """Best-effort default outputs root.

    Order: ``$STRIDE_OUTPUTS_ROOT`` → the current directory if it already holds
    ``outputs_s7/`` → the parent directory (the common case where the analysis
    repo sits next to this project).
    """
    env = os.environ.get("STRIDE_OUTPUTS_ROOT")
    if env:
        return Path(env)
    here = Path.cwd()
    if (here / "outputs_s7").is_dir():
        return here
    if (here.parent / "outputs_s7").is_dir():
        return here.parent
    return here


def resolve_paths(
    root: str | os.PathLike[str] | None = None,
    output: str | os.PathLike[str] | None = None,
) -> Paths:
    """Resolve the outputs root and the figure output directory."""
    root_path = Path(root) if root is not None else default_root()
    out_path = (
        Path(output)
        if output is not None
        else Path(__file__).resolve().parent / "output"
    )
    out_path.mkdir(parents=True, exist_ok=True)
    return Paths(root=root_path, output=out_path)
