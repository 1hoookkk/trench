"""Batch pill sampler — generate N admitted pill candidates into _queue/.

Draws from a scope envelope (via reroll.py), compiles each candidate
(compile_pill.py), and runs pill_gate.py admission. Only admitted pills land
in the queue. Candidates that fail gating are discarded.

Bias (enabled by default when verdicts.csv exists):
  Reads verdicts.csv, computes a weighted centroid shift of the scope
  envelope, and splits the batch into an exploit fraction (biased envelope)
  and an explore fraction (unshifted envelope). The shift formula and weights
  are locked in tools/workbench/bias.py.

Usage (from repo root):
    python tools/workbench/sampler.py --scope hf_resonance --n 20
    python tools/workbench/sampler.py --n 50 --seed 42
    python tools/workbench/sampler.py --fresh --n 20      # wipe queue, fresh batch
    python tools/workbench/sampler.py --no-bias --n 20    # disable centroid shift

Queue is persistent across sessions. Run verdict.py with --resume to skip
pills already logged.
"""
from __future__ import annotations

import argparse
import json
import math
import random
import sys
from pathlib import Path

_repo = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_repo / "tools"))
sys.path.insert(0, str(_repo / "tools" / "cube_authoring"))

from reroll import sample_pill, SCOPE_ENVELOPES          # noqa: E402
from compile_pill import compile_pill                     # noqa: E402
from pill_gate import run_pill_gate                       # noqa: E402
from bias import (                                        # noqa: E402
    extract_bias_params, pill_fingerprint,
    load_verdicts, shift_envelope, EXPLORE_FRAC,
)

WORKBENCH_DIR = Path(__file__).parent
QUEUE_DIR = WORKBENCH_DIR / "_queue"
VERDICTS_CSV = WORKBENCH_DIR / "verdicts.csv"


def _fail_reasons(verdict: dict) -> str:
    reasons: list[str] = []
    for p in verdict.get("per_probe", []):
        if p.get("status") == "failed":
            for e in p.get("edges", []):
                reasons.extend(e.get("fails", []))
    return ", ".join(sorted(set(reasons))) or "unknown"


def _load_queue_fingerprints(queue_dir: Path) -> set[str]:
    fps: set[str] = set()
    for fpath in queue_dir.glob("*.json"):
        try:
            pill = json.loads(fpath.read_text(encoding="utf-8"))
            bp = pill.get("provenance", {}).get("bias_params")
            if bp:
                fps.add(pill_fingerprint(bp))
        except Exception:
            pass
    return fps


def _print_bias_report(shift_report: dict, explore_n: int, exploit_n: int) -> None:
    n = shift_report["n_verdicts"]
    total = explore_n + exploit_n
    explore_pct = explore_n / total * 100 if total else 0.0
    print(f"bias: {n} verdicts  explore={explore_n}/{total} ({explore_pct:.0f}%)")
    for corner_key, cr in shift_report.get("corners", {}).items():
        print(
            f"  {corner_key}: "
            f"freq {cr['freq_center_before']}→{cr['freq_center_after']} "
            f"({cr['freq_shift_pct']:+.1f}%)  "
            f"gain {cr['gain_center_before']}→{cr['gain_center_after']} "
            f"({cr['gain_shift_pct']:+.1f}%)"
        )
    print()


def _generate_batch(
    scope_name: str,
    scope_env: dict,
    target: int,
    label: str,
    rng: random.Random,
    out_dir: Path,
    fingerprints: set[str],
    seq_start: int,
    admitted_offset: int,
    total_target: int,
) -> tuple[int, int, int]:
    """Generate `target` admitted pills from `scope_env`.

    Returns (admitted, attempts, next_seq).
    Mutates `fingerprints` in place to suppress duplicates within this run.
    """
    admitted = attempts = 0
    seq = seq_start
    dupes = 0

    while admitted < target:
        attempts += 1
        pill_id = f"{scope_name}.workbench.{label}.{attempts:04d}"
        authored = sample_pill(scope_name, rng, pill_id=pill_id, name=pill_id,
                               scope_env=scope_env)

        bp = extract_bias_params(authored)
        fp = pill_fingerprint(bp)
        if fp in fingerprints:
            dupes += 1
            continue
        fingerprints.add(fp)

        compiled = compile_pill(authored)
        compiled.setdefault("provenance", {})["bias_params"] = bp

        verdict = run_pill_gate(compiled)
        if not verdict["admitted"]:
            reasons = _fail_reasons(verdict)
            print(f"  [    ] rejected  {label}  attempt {attempts:5d}  [{reasons}]")
            continue

        admitted += 1
        global_n = admitted_offset + admitted
        fname = out_dir / f"{scope_name}_{seq:04d}.json"
        fname.write_text(json.dumps(compiled, indent=2) + "\n", encoding="utf-8")
        seq += 1
        print(f"  [{global_n:4d}/{total_target}] admitted  {label}  attempt {attempts:5d}"
              + (f"  (dupe skip×{dupes})" if dupes else ""))
        dupes = 0

    return admitted, attempts, seq


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
    parser.add_argument("--no-bias", action="store_true",
                        help="Disable centroid shift even if verdicts.csv exists.")
    parser.add_argument("--fresh", "--wipe-queue", action="store_true",
                        help="Wipe _queue/ before generating. Starts seq from 0.")
    args = parser.parse_args(argv)

    out_dir: Path = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.fresh:
        wiped = 0
        for f in sorted(out_dir.glob(f"{args.scope}_*.json")):
            f.unlink()
            wiped += 1
        print(f"--fresh: wiped {wiped} existing pills from {out_dir}")

    existing = sorted(out_dir.glob(f"{args.scope}_*.json"))
    seq = len(existing)
    fingerprints = _load_queue_fingerprints(out_dir)

    rng = random.Random(args.seed)

    # --- Bias ---
    original_env = SCOPE_ENVELOPES[args.scope]
    bias_active = not args.no_bias
    shifted_env = original_env
    shift_report: dict = {"n_verdicts": 0, "corners": {}}

    if bias_active and VERDICTS_CSV.exists():
        verdicts = load_verdicts(VERDICTS_CSV, args.scope)
        if verdicts:
            shifted_env, shift_report = shift_envelope(original_env, verdicts)

    has_bias = bool(shift_report.get("n_verdicts"))
    explore_n = math.ceil(args.n * EXPLORE_FRAC) if has_bias else 0
    exploit_n = args.n - explore_n

    print(f"scope={args.scope!r}  target={args.n} pills  seed={args.seed}")
    print(f"queue: {out_dir}  (seq starts at {seq}, {len(existing)} existing, "
          f"{len(fingerprints)} fingerprints)")
    if has_bias:
        _print_bias_report(shift_report, explore_n, exploit_n)
    else:
        print(f"bias: off (no verdicts)  explore=0/{args.n}")
        print()

    # --- Generate exploit batch (shifted/biased envelope) ---
    total_admitted = 0
    total_attempts = 0

    if exploit_n > 0:
        adm, att, seq = _generate_batch(
            args.scope, shifted_env, exploit_n,
            label="exploit", rng=rng, out_dir=out_dir,
            fingerprints=fingerprints, seq_start=seq,
            admitted_offset=0, total_target=args.n,
        )
        total_admitted += adm
        total_attempts += att

    # --- Generate explore batch (unshifted envelope) ---
    if explore_n > 0:
        adm, att, seq = _generate_batch(
            args.scope, original_env, explore_n,
            label="explore", rng=rng, out_dir=out_dir,
            fingerprints=fingerprints, seq_start=seq,
            admitted_offset=total_admitted, total_target=args.n,
        )
        total_admitted += adm
        total_attempts += att

    rate = total_admitted / total_attempts * 100 if total_attempts else 0.0
    realized_explore_pct = explore_n / args.n * 100 if args.n else 0.0
    print(f"\nadmitted {total_admitted}/{total_attempts} ({rate:.0f}% pass rate)  "
          f"exploration={realized_explore_pct:.0f}%")
    print(f"queue: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
