"""Shared helpers: data loading, deterministic saving, and table export.

Loaders read the pipeline's parquet outputs (never raw STRIDE). Savers write every
figure as SVG + PDF + PNG(600 dpi) with timestamps stripped so output is
byte-reproducible. Table exporters emit booktabs LaTeX, HTML, and Markdown.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.figure import Figure

from figure_config import CATALYTIC_DOMAINS, DPI, FORMATS, Paths

# --------------------------------------------------------------------------
# Loading
# --------------------------------------------------------------------------
def load_s7_figure(paths: Paths, stem: str) -> pd.DataFrame:
    """Load an S7 prepared-figure table, e.g. ``F1_reproducibility_landscape``."""
    path = paths.s7(f"{stem}.parquet")
    if not path.is_file():
        raise FileNotFoundError(
            f"missing S7 figure input: {path}\n"
            "Run the STRIDE pipeline through S7 first (see README)."
        )
    return pd.read_parquet(path)


def load_table(paths: Paths, stage: str, stem: str) -> pd.DataFrame:
    """Load any stage table, e.g. ``load_table(paths, 's4', 'significance_screen')``."""
    path = paths.table(stage, f"{stem}.parquet")
    if not path.is_file():
        raise FileNotFoundError(f"missing input table: {path}")
    return pd.read_parquet(path)


def is_catalytic(domain: pd.Series) -> pd.Series:
    """Boolean mask: rows whose domain is catalytic machinery."""
    return domain.astype(str).isin(CATALYTIC_DOMAINS)


def parse_canon_label(label: str) -> tuple[str, int]:
    """Split a canonical label ``"{chain}:{resid}"`` into ``(chain, resid)``.

    Falls back to ``(label, 0)`` for any label without the ``:`` separator so
    ordering never raises on unexpected input.
    """
    s = str(label)
    if ":" not in s:
        return (s, 0)
    chain, _, resid = s.rpartition(":")
    try:
        return (chain, int(resid))
    except ValueError:
        return (chain, 0)


def residue_order(labels: pd.Series) -> list[str]:
    """Unique canonical labels ordered by (chain, residue number)."""
    uniq = list(dict.fromkeys(labels.astype(str)))
    return sorted(uniq, key=parse_canon_label)


# --------------------------------------------------------------------------
# Deterministic multi-format saving
# --------------------------------------------------------------------------
_METADATA = {
    "svg": {"Date": None},
    "pdf": {"Creator": "paper_figures", "Producer": "paper_figures",
            "CreationDate": None},
    "png": {"Software": None},
}


def save_figure(fig: Figure, paths: Paths, stem: str) -> list[Path]:
    """Write ``fig`` to ``output/<stem>.{svg,pdf,png}`` deterministically."""
    written: list[Path] = []
    for ext in FORMATS:
        out = paths.output / f"{stem}.{ext}"
        fig.savefig(
            out,
            format=ext,
            dpi=DPI if ext == "png" else None,
            metadata=_METADATA.get(ext),
        )
        written.append(out)
    plt.close(fig)
    return written


def panel_letter(ax: plt.Axes, letter: str, *, dx: float = -0.02, dy: float = 1.04) -> None:
    """Place a bold panel letter (A, B, …) at the top-left of an axes."""
    ax.text(
        dx,
        dy,
        letter,
        transform=ax.transAxes,
        fontsize=10,
        fontweight="bold",
        va="bottom",
        ha="right",
    )


def order_categorical(values: pd.Series, order: tuple[str, ...]) -> list[str]:
    """Return the members of ``order`` present in ``values`` (stable), then extras."""
    present = set(values.astype(str))
    ordered = [v for v in order if v in present]
    extras = sorted(present - set(order))
    return ordered + extras


# --------------------------------------------------------------------------
# Manuscript table export (LaTeX booktabs / HTML / Markdown)
# --------------------------------------------------------------------------
def _fmt_cell(v: object) -> str:
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, float):
        if v != v:  # NaN
            return "--"
        return f"{v:.3g}"
    return str(v)


def _escape_latex(s: str) -> str:
    repl = {
        "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#", "_": r"\_",
        "{": r"\{", "}": r"\}", "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    for a, b in repl.items():
        s = s.replace(a, b)
    return s


def export_table(
    df: pd.DataFrame,
    paths: Paths,
    stem: str,
    *,
    caption: str,
    label: str,
) -> list[Path]:
    """Export a table as ``.tex`` (booktabs), ``.html`` and ``.md``."""
    written: list[Path] = []
    cols = list(df.columns)
    rows = [[_fmt_cell(v) for v in row] for row in df.itertuples(index=False, name=None)]

    # ---- LaTeX (booktabs) ------------------------------------------------
    colspec = "l" + "r" * (len(cols) - 1) if cols else "l"
    tex = [
        r"\begin{table}[t]",
        r"  \centering",
        rf"  \caption{{{_escape_latex(caption)}}}",
        rf"  \label{{{label}}}",
        rf"  \begin{{tabular}}{{{colspec}}}",
        r"    \toprule",
        "    " + " & ".join(_escape_latex(c.replace("_", " ")) for c in cols) + r" \\",
        r"    \midrule",
    ]
    for r in rows:
        tex.append("    " + " & ".join(_escape_latex(c) for c in r) + r" \\")
    tex += [r"    \bottomrule", r"  \end{tabular}", r"\end{table}", ""]
    tex_path = paths.output / f"{stem}.tex"
    tex_path.write_text("\n".join(tex))
    written.append(tex_path)

    # ---- HTML ------------------------------------------------------------
    html = [
        f'<table class="stride-table"><caption>{caption}</caption>',
        "<thead><tr>" + "".join(f"<th>{c.replace('_', ' ')}</th>" for c in cols)
        + "</tr></thead>",
        "<tbody>",
    ]
    for r in rows:
        html.append("<tr>" + "".join(f"<td>{c}</td>" for c in r) + "</tr>")
    html += ["</tbody></table>", ""]
    html_path = paths.output / f"{stem}.html"
    html_path.write_text("\n".join(html))
    written.append(html_path)

    # ---- Markdown --------------------------------------------------------
    md = [
        f"**{caption}**",
        "",
        "| " + " | ".join(c.replace("_", " ") for c in cols) + " |",
        "| " + " | ".join("---" for _ in cols) + " |",
    ]
    for r in rows:
        md.append("| " + " | ".join(r) + " |")
    md.append("")
    md_path = paths.output / f"{stem}.md"
    md_path.write_text("\n".join(md))
    written.append(md_path)

    return written
