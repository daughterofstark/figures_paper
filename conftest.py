"""Pytest bootstrap: put the project root (and tests/) on sys.path.

The figure modules use flat imports (``import styles``), so the project directory
must be importable when the tests run.
"""
from __future__ import annotations

import os
import sys

_HERE = os.path.dirname(os.path.abspath(__file__))
for _p in (_HERE, os.path.join(_HERE, "tests")):
    if _p not in sys.path:
        sys.path.insert(0, _p)
