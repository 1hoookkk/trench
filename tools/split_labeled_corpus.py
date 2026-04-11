"""Create deterministic, class-balanced design/holdout corpus splits.

Uses label classification from tools.promote_better_than_emu.classify_name.
"""
from __future__ import annotations

import argparse
import random
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

import tools.promote_better_than_emu as scorer


def _clear_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for wav in path.glob("*.wav"):
        wav.unlink()


def main() -> None:
    parser = argparse.ArgumentParser(description="Split corpus into class-balanced design/holdout sets.")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--holdout", type=Path, required=True)
    parser.add_argument("--holdout-ratio", type=float, default=0.2)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--max-per-class",
        type=int,
        default=None,
        help="Optional cap applied per label before design/holdout split.",
    )
    args = parser.parse_args()

    if not args.source.exists():
        raise FileNotFoundError(f"Source corpus not found: {args.source}")
    if args.holdout_ratio <= 0.0 or args.holdout_ratio >= 1.0:
        raise ValueError("--holdout-ratio must be between 0 and 1.")
    if args.max_per_class is not None and args.max_per_class <= 0:
        raise ValueError("--max-per-class must be positive when provided.")

    _clear_dir(args.design)
    _clear_dir(args.holdout)

    grouped: dict[str, list[Path]] = {"bass": [], "vocal": [], "drum": [], "mix": []}
    unknown: list[Path] = []
    for wav in args.source.rglob("*.wav"):
        label = scorer.classify_name(wav)
        if label is None:
            unknown.append(wav)
        else:
            grouped[label].append(wav)

    rng = random.Random(args.seed)
    source_count = {k: len(v) for k, v in grouped.items()}
    capped_count = {k: 0 for k in grouped}
    design_count = {k: 0 for k in grouped}
    holdout_count = {k: 0 for k in grouped}

    for label, files in grouped.items():
        files = files[:]
        rng.shuffle(files)
        if args.max_per_class is not None:
            files = files[:args.max_per_class]
        capped_count[label] = len(files)
        n_holdout = max(1, int(round(len(files) * args.holdout_ratio))) if files else 0
        holdout_files = files[:n_holdout]
        design_files = files[n_holdout:]
        for i, src in enumerate(holdout_files, start=1):
            dst = args.holdout / f"{label}_{i:04d}_{src.name}"
            shutil.copy2(src, dst)
            holdout_count[label] += 1
        for i, src in enumerate(design_files, start=1):
            dst = args.design / f"{label}_{i:04d}_{src.name}"
            shutil.copy2(src, dst)
            design_count[label] += 1

    print(f"design={args.design}")
    print(f"holdout={args.holdout}")
    print(f"source_count={source_count}")
    print(f"capped_count={capped_count}")
    print(f"design_count={design_count}")
    print(f"holdout_count={holdout_count}")
    print(f"unknown_ignored={len(unknown)}")


if __name__ == "__main__":
    main()
