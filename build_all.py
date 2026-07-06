"""Regenerate every publication figure and manuscript table.

Usage::

    python build_all.py [--root <outputs-root>] [--output-dir <dir>]

``--root`` is the directory that contains ``outputs_s2/`` … ``outputs_s7/`` (the
completed STRIDE downstream run). If omitted it is auto-detected (see
``figure_config.default_root``). All artifacts are written under ``--output-dir``
(default ``paper_figures/output``) as SVG + PDF + PNG(600 dpi) for figures and
LaTeX + HTML + Markdown for tables.

This project never reads raw STRIDE data, never reruns S0–S7, and computes no new
statistics — it only visualises tables the pipeline already produced.
"""
from __future__ import annotations

import argparse
import importlib
import time
from pathlib import Path

from figure_config import Paths, resolve_paths

FIGURE_MODULES = [f"figure{i}" for i in range(1, 7)]
SUPPLEMENTARY_MODULES = [f"supplementary_figure{i}" for i in range(1, 3)]
TABLE_MODULE = "tables"


def build_all(paths: Paths) -> dict[str, list[Path]]:
    results: dict[str, list[Path]] = {}
    for name in [*FIGURE_MODULES, *SUPPLEMENTARY_MODULES, TABLE_MODULE]:
        mod = importlib.import_module(name)
        t0 = time.perf_counter()
        written = mod.build(paths)
        dt = time.perf_counter() - t0
        results[name] = written
        print(f"  {name:<22s} → {len(written):2d} files  ({dt:4.1f}s)")
    return results


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Build all paper figures + tables.")
    parser.add_argument("--root", default=None,
                        help="outputs root containing outputs_s2/ … outputs_s7/")
    parser.add_argument("--output-dir", default=None,
                        help="where to write figures/tables (default: ./output)")
    args = parser.parse_args(argv)

    paths = resolve_paths(root=args.root, output=args.output_dir)
    print(f"outputs root : {paths.root}")
    print(f"output dir   : {paths.output}")
    print("building:")
    results = build_all(paths)
    total = sum(len(v) for v in results.values())
    print(f"done: {total} files across {len(results)} modules → {paths.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
