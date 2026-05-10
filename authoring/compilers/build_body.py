"""Build one or more DSL body files: author .py -> raw JSON -> compiled JSON -> plots.

Usage:
  python tools/build_body.py authoring/bodies/trap_vowel_slam.py [...more.py]
"""
from __future__ import annotations

import argparse
import importlib.util
import subprocess
import sys
from pathlib import Path

_TOOLS = Path(__file__).resolve().parent
if str(_TOOLS) not in sys.path:
    sys.path.insert(0, str(_TOOLS))
_ROOT = _TOOLS.parent
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from body_dsl import Body, build  # noqa: E402


def load_body(path: Path) -> Body:
    spec = importlib.util.spec_from_file_location(path.stem, path)
    if spec is None or spec.loader is None:
        raise SystemExit(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    if not hasattr(module, "BODY"):
        raise SystemExit(f"{path} has no module-level BODY")
    body = module.BODY
    if not (hasattr(body, "name") and hasattr(body, "slots") and hasattr(body, "boost")):
        raise SystemExit(f"{path}.BODY missing expected Body fields (name/slots/boost)")
    return body


def main() -> int:
    ap = argparse.ArgumentParser(description="Build DSL body .py files into raw+compiled+plots.")
    ap.add_argument("bodies", nargs="+", type=Path, help="body .py files")
    ap.add_argument("--out-root", type=Path, default=_ROOT / "cartridges" / "engine" / "_source" / "body_dsl")
    ap.add_argument("--plot-dir", type=Path, default=_ROOT / "parity_plots")
    args = ap.parse_args()

    for body_path in args.bodies:
        body = load_body(body_path)
        raw_path, compiled_path = build(body, args.out_root)
        print(f"[build] {body.name}")
        print(f"  raw:      {raw_path}")
        print(f"  compiled: {compiled_path}")

        pair_png = args.plot_dir / f"{body.name}_stage_pairs.png"
        accum_png = args.plot_dir / f"{body.name}_stage_accumulation.png"
        subprocess.run([sys.executable, str(_TOOLS / "plot_stage_pairs.py"),
                        str(compiled_path), "--out", str(pair_png)], check=True)
        subprocess.run([sys.executable, str(_TOOLS / "plot_stage_accumulation.py"),
                        str(compiled_path), "--out", str(accum_png)], check=True)
        print(f"  pairs:    {pair_png}")
        print(f"  accum:    {accum_png}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
