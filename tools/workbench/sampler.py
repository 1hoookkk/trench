"""Batch pill sampler — generate N admitted pill candidates into _queue/.

Draws from a scope envelope (via reroll.py), compiles each candidate
(compile_pill.py), and runs pill_gate.py admission. Only admitted pills land
in the queue. Candidates that fail gating are discarded.

Usage (from repo root):
    python tools/workbench/sampler.py --scope hf_resonance --n 20
    python tools/workbench/sampler.py --n 50 --seed 42

Each admitted pill is written as:
    tools/workbench/_queue/<scope>_<seq04d>.json

Run verdict.py after to audition and log keystrokes.

Scope caveat: SCOPE_ENVELOPES in reroll.py has one entry today (hf_resonance).
Each new body scope needs an envelope authored there before this sampler can
generate candidates for it.
"""
from __future__ import annotations

import argparse
import json
import random
import sys
from pathlib import Path

_repo = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_repo / "tools"))
sys.path.insert(0, str(_repo / "tools" / "cube_authoring"))

from reroll import sample_pill, SCOPE_ENVELOPES          # noqa: E402
from compile_pill import compile_pill                     # noqa: E402
from pill_gate import run_pill_gate                       # noqa: E402

QUEUE_DIR = Path(__file__).parent / "_queue"


def _fail_reasons(verdict: dict) -> str:
    reasons: list[str] = []
    for p in verdict.get("per_probe", []):
        if p.get("status") == "failed":
            for e in p.get("edges", []):
                reasons.extend(e.get("fails", []))
    return ", ".join(sorted(set(reasons))) or "unknown"


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Generate N admitted pill candidates into _queue/.",
    )
    parser.add_argument("--scope", default="hf_resonance",
                        choices=list(SCOPE_ENVELOPES),
                        help="Scope name (default: hf_resonance).")
    parser.add_argument("--n", type=int, default=20,
                        help="Number of admitted pills to generate (default: 20).")
    parser.add_argument("--seed", type=int, default=None,
                        help="RNG seed for reproducibility (default: random).")
    parser.add_argument("--out-dir", type=Path, default=QUEUE_DIR,
                        help="Output directory (default: tools/workbench/_queue/).")
    args = parser.parse_args(argv)

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    rng = random.Random(args.seed)
    admitted = 0
    attempts = 0

    existing = sorted(out_dir.glob(f"{args.scope}_*.json"))
    seq = len(existing)

    print(f"scope={args.scope!r}  target={args.n} pills  seed={args.seed}")
    print(f"queue: {out_dir}  (seq starts at {seq}, {len(existing)} existing)")
    print()

    while admitted < args.n:
        attempts += 1
        pill_id = f"{args.scope}.workbench.{attempts:04d}"
        authored = sample_pill(args.scope, rng, pill_id=pill_id, name=pill_id)
        compiled = compile_pill(authored)
        verdict = run_pill_gate(compiled)

        if verdict["admitted"]:
            admitted += 1
            fname = out_dir / f"{args.scope}_{seq:04d}.json"
            fname.write_text(json.dumps(compiled, indent=2) + "\n", encoding="utf-8")
            seq += 1
            print(f"  [{admitted:4d}/{args.n}] admitted    attempt {attempts:5d}")
        else:
            reasons = _fail_reasons(verdict)
            print(f"  [    ] rejected    attempt {attempts:5d}  [{reasons}]")

    rate = admitted / attempts * 100 if attempts else 0.0
    print(f"\nadmitted {admitted}/{attempts} ({rate:.0f}% pass rate)")
    print(f"queue: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
