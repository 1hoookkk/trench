#!/usr/bin/env python3
"""Root launcher for the EOS/TRENCH encoded-surface compiler harness."""

from __future__ import annotations

import runpy
from pathlib import Path


if __name__ == "__main__":
    target = Path(__file__).resolve().parent / "authoring" / "compilers" / "trench_eos_pipeline_compiler.py"
    runpy.run_path(str(target), run_name="__main__")
