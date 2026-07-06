"""Publication styling: typography, colour, and determinism for paper_figures.

Centralises every visual decision so all figures share one identity:

* a clean sans-serif type scale with no top/right spines and no grid clutter,
* the Okabe–Ito colour-blind-safe qualitative palette and perceptually-uniform
  sequential maps (``viridis`` / ``cividis``) for continuous fields,
* deterministic output (fixed SVG hash-salt; no timestamps embedded on save).

Nothing here reads data or performs analysis — it only configures matplotlib.
"""
from __future__ import annotations

import matplotlib as mpl
import matplotlib.pyplot as plt

# --------------------------------------------------------------------------
# Colour-blind-safe qualitative palette (Okabe & Ito, 2008)
# --------------------------------------------------------------------------
OKABE_ITO = {
    "blue": "#0072B2",
    "orange": "#E69F00",
    "green": "#009E73",
    "vermillion": "#D55E00",
    "sky": "#56B4E9",
    "reddish_purple": "#CC79A7",
    "yellow": "#F0E442",
    "black": "#000000",
    "grey": "#999999",
}

# --------------------------------------------------------------------------
# SEMANTIC COLOUR SYSTEM  —  one colour = one meaning, used identically in every
# figure. Do not reuse a colour for an unrelated concept.
#
#   catalytic machinery ....... vermillion  (rings / accents ONLY)
#   serotype identity ......... blue / sky / green / purple  (DENV1..4)
#   reproducibility ρ (cont.).. viridis sequential
#   cross-serotype conservation cividis / blue sequential
#   FDR significance .......... darkness: significant = ink, else light grey
#   variance τ² (limiting) .... orange ;  σ̄² (sampling) = neutral grey
#   correlation (S-matrix) .... cividis diverging-through-neutral
# --------------------------------------------------------------------------

#: Catalytic machinery — reserved; never used for anything else.
CATALYTIC_ACCENT = OKABE_ITO["vermillion"]

#: Serotype identity (DENV1–4). Deliberately excludes orange/vermillion so the
#: catalytic and variance semantics stay unambiguous.
SEROTYPE_COLORS = {
    "DENV1": OKABE_ITO["blue"],
    "DENV2": OKABE_ITO["sky"],
    "DENV3": OKABE_ITO["green"],
    "DENV4": OKABE_ITO["reddish_purple"],
}

#: Cross-serotype conservation class (ordinal, dark = more conserved).
CONSERVATION_COLORS = {
    "reproducible_all": "#08519c",
    "reproducible_majority": "#3182bd",
    "reproducible_some": "#9ecae1",
    "reproducible_none": "#deebf7",
}

#: Variance components. τ² (between-replicate) is the reproducibility-limiting
#: term and gets the salient warm colour; σ̄² (sampling) is neutral.
VARIANCE_COLORS = {
    "tau2": OKABE_ITO["orange"],
    "sigma2": "#BDBDBD",
}

#: FDR significance encoded by darkness (no hue clash with catalytic/variance).
SIG_COLOR = "#222222"
NONSIG_COLOR = "#BBBBBB"

#: Perceptually-uniform sequential / correlation maps (colour-blind friendly).
SEQ_CMAP = "viridis"
SEQ_CMAP_ALT = "cividis"
CORR_CMAP = "cividis"

# --------------------------------------------------------------------------
# Type scale (points) — consistent across every panel
# --------------------------------------------------------------------------
FS_TITLE = 8
FS_LABEL = 7.5
FS_TICK = 6.5
FS_ANNOT = 6
FS_PANEL = 10  # panel letters (A, B, C …)


def set_determinism() -> None:
    """Make figure bytes reproducible run-to-run (fixed salt, no timestamps)."""
    mpl.rcParams["svg.hashsalt"] = "stride-paper-figures"


def apply_style() -> None:
    """Apply the global publication rcParams. Call once before plotting."""
    set_determinism()
    plt.rcParams.update(
        {
            # ---- typography -------------------------------------------------
            "font.family": "sans-serif",
            "font.sans-serif": [
                "Arial",
                "Helvetica",
                "DejaVu Sans",
            ],
            "font.size": FS_LABEL,
            "axes.titlesize": FS_TITLE,
            "axes.labelsize": FS_LABEL,
            "xtick.labelsize": FS_TICK,
            "ytick.labelsize": FS_TICK,
            "legend.fontsize": FS_ANNOT,
            "figure.titlesize": FS_TITLE,
            # keep text as editable text in vector output (not paths)
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            # ---- spare, journal-style axes ---------------------------------
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.linewidth": 0.6,
            "axes.grid": False,
            "axes.axisbelow": True,
            "xtick.direction": "out",
            "ytick.direction": "out",
            "xtick.major.width": 0.6,
            "ytick.major.width": 0.6,
            "xtick.major.size": 2.5,
            "ytick.major.size": 2.5,
            "lines.linewidth": 1.0,
            "lines.markersize": 4,
            "legend.frameon": False,
            "legend.handlelength": 1.2,
            "legend.columnspacing": 1.0,
            "legend.labelspacing": 0.3,
            # ---- save defaults ---------------------------------------------
            "figure.dpi": 150,
            "savefig.dpi": 600,
            "savefig.bbox": "tight",
            "savefig.pad_inches": 0.02,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def serotype_color(serotype: str) -> str:
    """Colour for a serotype, falling back to grey for unknown labels."""
    return SEROTYPE_COLORS.get(str(serotype), OKABE_ITO["grey"])
