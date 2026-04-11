"""Build shipping name targets from P2k-trained role vocabulary.

Targets are keyed by shipping names and sourced from mapped P2k bodies.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

DEFAULT_MAP = {
    "Speaker Knockerz": 6,
    "Aluminum Siding": 26,
    "Small Talk Ah-Ee": 13,
    "Cul-De-Sac": 10,
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Build shipping role targets from P2k vocabulary.")
    parser.add_argument(
        "--vocab",
        type=Path,
        default=ROOT / "datasets" / "role_vocab" / "p2k_role_vocabulary_v1.json",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "datasets" / "role_vocab" / "shipping_role_targets_v1.json",
    )
    args = parser.parse_args()

    if not args.vocab.exists():
        raise FileNotFoundError(f"Vocabulary file not found: {args.vocab}")

    vocab = json.loads(args.vocab.read_text(encoding="utf-8"))
    by_key = {row["key"]: row for row in vocab.get("bodies", [])}
    targets = {}

    for shipping_name, idx in DEFAULT_MAP.items():
        key = f"P2k_{idx:03d}"
        if key not in by_key:
            raise KeyError(f"Missing {key} in {args.vocab}")
        targets[shipping_name] = {
            "source_p2k": key,
            "source_name": by_key[key]["name"],
            "signature": by_key[key]["signature"],
        }

    payload = {
        "version": "shipping-role-targets-v1",
        "vocab_version": vocab.get("version"),
        "targets": targets,
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(args.out)
    print(f"targets={list(targets.keys())}")


if __name__ == "__main__":
    main()
