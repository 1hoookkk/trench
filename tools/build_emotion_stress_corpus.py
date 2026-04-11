"""Build an emotion-stress vocal corpus and merge with non-vocal classes.

Emotion proxy buckets (from VocalSet filename tokens):
- intense: belt, fry
- fragile: breathy, inhaled, pp
- unstable: vibrato, trill, trillo, messa
- neutral: straight, spoken
"""
from __future__ import annotations

import argparse
import json
import random
import re
import shutil
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

import tools.promote_better_than_emu as scorer


BUCKET_TOKENS: dict[str, tuple[str, ...]] = {
    "intense": ("belt", "fry"),
    "fragile": ("breathy", "inhaled", "pp"),
    "unstable": ("vibrato", "trill", "trillo", "messa"),
    "neutral": ("straight", "spoken"),
}


def _clear_wavs(path: Path) -> None:
    path.mkdir(parents=True, exist_ok=True)
    for wav in path.glob("*.wav"):
        wav.unlink()


def _tokenize(path: Path) -> set[str]:
    return set(t for t in re.split(r"[^a-z0-9]+", path.stem.lower()) if t)


def _collect_bucket_candidates(vocal_source: Path) -> dict[str, list[Path]]:
    candidates: dict[str, list[Path]] = {k: [] for k in BUCKET_TOKENS}
    for wav in vocal_source.rglob("*.wav"):
        tokens = _tokenize(wav)
        for bucket, bucket_tokens in BUCKET_TOKENS.items():
            if any(tok in tokens for tok in bucket_tokens):
                candidates[bucket].append(wav)
    return candidates


def _pick_bucket_files(
    candidates: dict[str, list[Path]],
    *,
    per_bucket: int,
    seed: int,
) -> dict[str, list[Path]]:
    rng = random.Random(seed)
    selected: dict[str, list[Path]] = {k: [] for k in BUCKET_TOKENS}
    used: set[Path] = set()

    order = sorted(BUCKET_TOKENS.keys(), key=lambda b: len(candidates[b]))
    for bucket in order:
        pool = candidates[bucket][:]
        rng.shuffle(pool)
        chosen: list[Path] = []
        for path in pool:
            if path in used:
                continue
            chosen.append(path)
            used.add(path)
            if len(chosen) >= per_bucket:
                break
        if len(chosen) < per_bucket:
            raise RuntimeError(
                f"Not enough unique files for bucket '{bucket}'. "
                f"Needed {per_bucket}, got {len(chosen)}."
            )
        selected[bucket] = chosen
    return selected


def _copy_selected_vocals(selected: dict[str, list[Path]], out_vocals: Path) -> int:
    _clear_wavs(out_vocals)
    total = 0
    for bucket in BUCKET_TOKENS:
        files = selected[bucket]
        for i, src in enumerate(files, start=1):
            dst = out_vocals / f"vocal_emotion_{bucket}_{i:04d}_{src.name}"
            shutil.copy2(src, dst)
            total += 1
    return total


def _merge_with_non_vocals(base_source: Path, out_vocals: Path, out_corpus: Path) -> dict[str, int]:
    _clear_wavs(out_corpus)

    for wav in base_source.rglob("*.wav"):
        label = scorer.classify_name(wav)
        if label is None or label == "vocal":
            continue
        dst = out_corpus / wav.name
        shutil.copy2(wav, dst)

    for wav in out_vocals.glob("*.wav"):
        dst = out_corpus / wav.name
        shutil.copy2(wav, dst)

    return scorer.count_corpus_labels(out_corpus)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build emotion-stress vocal corpus.")
    parser.add_argument(
        "--vocal-source",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "vocals_vocalset",
    )
    parser.add_argument(
        "--base-source",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "user_kits_real",
    )
    parser.add_argument(
        "--out-vocals",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "vocals_emotion_stress",
    )
    parser.add_argument(
        "--out-corpus",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "user_kits_plus_emotion_vocals",
    )
    parser.add_argument("--per-bucket", type=int, default=120)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    if args.per_bucket <= 0:
        raise ValueError("--per-bucket must be positive.")
    if not args.vocal_source.exists():
        raise FileNotFoundError(f"Missing vocal source: {args.vocal_source}")
    if not args.base_source.exists():
        raise FileNotFoundError(f"Missing base source: {args.base_source}")

    candidates = _collect_bucket_candidates(args.vocal_source)
    selected = _pick_bucket_files(candidates, per_bucket=args.per_bucket, seed=args.seed)
    total_vocal = _copy_selected_vocals(selected, args.out_vocals)
    merged_counts = _merge_with_non_vocals(args.base_source, args.out_vocals, args.out_corpus)

    manifest = {
        "vocal_source": str(args.vocal_source),
        "base_source": str(args.base_source),
        "out_vocals": str(args.out_vocals),
        "out_corpus": str(args.out_corpus),
        "seed": args.seed,
        "per_bucket": args.per_bucket,
        "candidate_counts": {k: len(v) for k, v in candidates.items()},
        "selected_counts": {k: len(v) for k, v in selected.items()},
        "merged_label_counts": merged_counts,
        "bucket_tokens": BUCKET_TOKENS,
    }
    manifest_path = args.out_vocals / "emotion_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    print(f"out_vocals={args.out_vocals}")
    print(f"out_corpus={args.out_corpus}")
    print(f"total_vocal={total_vocal}")
    print(f"candidate_counts={manifest['candidate_counts']}")
    print(f"selected_counts={manifest['selected_counts']}")
    print(f"merged_label_counts={merged_counts}")
    print(f"manifest={manifest_path}")


if __name__ == "__main__":
    main()
