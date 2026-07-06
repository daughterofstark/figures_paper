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


def chain_bands(chains: list[str]) -> list[tuple[str, int, int]]:
    """Contiguous (chain, start_idx, end_idx_exclusive) spans over an ordered list.

    Used to shade chain regions along an ordinal residue axis.
    """
    bands: list[tuple[str, int, int]] = []
    if not chains:
        return bands
    start = 0
    for i in range(1, len(chains) + 1):
        if i == len(chains) or chains[i] != chains[start]:
            bands.append((chains[start], start, i))
            start = i
    return bands


def lowess(
    x: "np.ndarray", y: "np.ndarray", *, frac: float = 0.6, n_out: int = 120
) -> "tuple[np.ndarray, np.ndarray] | None":
    """Deterministic LOWESS smoother (tricube weights, local linear fit).

    A visualization aid only — it draws a trend through points that already exist;
    it does not create or claim any new statistic about the data. Returns
    ``(xs, ys)`` or ``None`` when there are too few finite points.
    """
    import numpy as np

    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x, y = x[m], y[m]
    if len(x) < 4 or np.ptp(x) == 0:
        return None
    order = np.argsort(x, kind="stable")
    x, y = x[order], y[order]
    xs = np.linspace(x.min(), x.max(), n_out)
    r = max(2, int(np.ceil(frac * len(x))))
    ys = np.empty_like(xs)
    for i, xv in enumerate(xs):
        d = np.abs(x - xv)
        h = np.sort(d)[min(r - 1, len(d) - 1)] or 1.0
        w = np.clip(d / h, 0.0, 1.0)
        w = (1.0 - w**3) ** 3
        sw = w.sum()
        if sw <= 0:
            ys[i] = np.nan
            continue
        mx = (w * x).sum() / sw
        my = (w * y).sum() / sw
        bxx = (w * (x - mx) ** 2).sum()
        b = (w * (x - mx) * (y - my)).sum() / bxx if bxx > 0 else 0.0
        ys[i] = (my - b * mx) + b * xv
    return xs, ys


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
