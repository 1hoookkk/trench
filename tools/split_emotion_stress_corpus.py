"""Split corpus with strict vocal emotion-bucket stratification.

Expected vocal filenames include one of:
- emotion_intense
- emotion_fragile
- emotion_unstable
- emotion_neutral
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

VOCAL_BUCKETS = ("intense", "fragile", "unstable", "neutral")


def _clear_dir(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for wav in path.glob("*.wav"):
        wav.unlink()


def _vocal_bucket(path: Path) -> str | None:
    stem = path.stem.lower()
    for bucket in VOCAL_BUCKETS:
        if f"emotion_{bucket}" in stem:
            return bucket
    return None


def _alloc(total: int, n: int) -> list[int]:
    base = total // n
    rem = total % n
    return [base + (1 if i < rem else 0) for i in range(n)]


def _copy_many(files: list[Path], out_dir: Path, label: str, start_idx: int) -> int:
    idx = start_idx
    for src in files:
        dst = out_dir / f"{label}_{idx:04d}_{src.name}"
        shutil.copy2(src, dst)
        idx += 1
    return idx


def main() -> None:
    parser = argparse.ArgumentParser(description="Split corpus with vocal emotion-bucket stratification.")
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--design", type=Path, required=True)
    parser.add_argument("--holdout", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--max-per-class", type=int, default=70)
    parser.add_argument("--holdout-ratio", type=float, default=0.2)
    args = parser.parse_args()

    if not args.source.exists():
        raise FileNotFoundError(f"Source corpus not found: {args.source}")
    if args.max_per_class <= 0:
        raise ValueError("--max-per-class must be positive.")
    if args.holdout_ratio <= 0.0 or args.holdout_ratio >= 1.0:
        raise ValueError("--holdout-ratio must be between 0 and 1.")

    _clear_dir(args.design)
    _clear_dir(args.holdout)

    grouped: dict[str, list[Path]] = {"bass": [], "vocal": [], "drum": [], "mix": []}
    vocal_grouped: dict[str, list[Path]] = {b: [] for b in VOCAL_BUCKETS}
    unknown = 0
    for wav in args.source.rglob("*.wav"):
        label = scorer.classify_name(wav)
        if label is None:
            unknown += 1
            continue
        grouped[label].append(wav)
        if label == "vocal":
            b = _vocal_bucket(wav)
            if b is not None:
                vocal_grouped[b].append(wav)

    rng = random.Random(args.seed)
    source_count = {k: len(v) for k, v in grouped.items()}
    source_vocal_bucket_count = {k: len(v) for k, v in vocal_grouped.items()}

    design_count = {k: 0 for k in grouped}
    holdout_count = {k: 0 for k in grouped}

    # Non-vocal classes use normal capped split.
    for label in ("bass", "drum", "mix"):
        files = grouped[label][:]
        rng.shuffle(files)
        files = files[:args.max_per_class]
        n_holdout = max(1, int(round(len(files) * args.holdout_ratio))) if files else 0
        holdout_files = files[:n_holdout]
        design_files = files[n_holdout:]
        _copy_many(holdout_files, args.holdout, label, holdout_count[label] + 1)
        holdout_count[label] += len(holdout_files)
        _copy_many(design_files, args.design, label, design_count[label] + 1)
        design_count[label] += len(design_files)

    # Vocal class: strict per-bucket quotas.
    vocal_target_total = args.max_per_class
    vocal_holdout_total = max(1, int(round(vocal_target_total * args.holdout_ratio)))
    vocal_design_total = vocal_target_total - vocal_holdout_total

    holdout_bucket_quota = _alloc(vocal_holdout_total, len(VOCAL_BUCKETS))
    design_bucket_quota = _alloc(vocal_design_total, len(VOCAL_BUCKETS))

    design_vocal_bucket_count = {k: 0 for k in VOCAL_BUCKETS}
    holdout_vocal_bucket_count = {k: 0 for k in VOCAL_BUCKETS}

    for i, bucket in enumerate(VOCAL_BUCKETS):
        pool = vocal_grouped[bucket][:]
        rng.shuffle(pool)
        need_hold = holdout_bucket_quota[i]
        need_design = design_bucket_quota[i]
        need_total = need_hold + need_design
        if len(pool) < need_total:
            raise RuntimeError(
                f"Not enough files in vocal bucket '{bucket}'. Need {need_total}, have {len(pool)}."
            )
        hold_files = pool[:need_hold]
        design_files = pool[need_hold:need_hold + need_design]

        _copy_many(hold_files, args.holdout, "vocal", holdout_count["vocal"] + 1)
        holdout_count["vocal"] += len(hold_files)
        holdout_vocal_bucket_count[bucket] = len(hold_files)

        _copy_many(design_files, args.design, "vocal", design_count["vocal"] + 1)
        design_count["vocal"] += len(design_files)
        design_vocal_bucket_count[bucket] = len(design_files)

    print(f"design={args.design}")
    print(f"holdout={args.holdout}")
    print(f"source_count={source_count}")
    print(f"source_vocal_bucket_count={source_vocal_bucket_count}")
    print(f"design_count={design_count}")
    print(f"holdout_count={holdout_count}")
    print(f"design_vocal_bucket_count={design_vocal_bucket_count}")
    print(f"holdout_vocal_bucket_count={holdout_vocal_bucket_count}")
    print(f"unknown_ignored={unknown}")


if __name__ == "__main__":
    main()
