"""Rank shipping finalists by P2k-trained role vocabulary distance.

This does not use P2k as an audio-metric baseline. It uses P2k-derived
stage-role signatures as target vocabulary per shipping name.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from pyruntime.analysis import shipping_gate
from pyruntime.body import Body
from pyruntime.role_vocab import body_signature, role_distance


BODY_KEYS = {
    "speaker_knockerz": "Speaker Knockerz",
    "aluminum_siding": "Aluminum Siding",
    "small_talk": "Small Talk Ah-Ee",
    "cul_de_sac": "Cul-De-Sac",
}


def main() -> None:
    parser = argparse.ArgumentParser(description="Promote by role vocabulary distance.")
    parser.add_argument(
        "--manifest",
        type=Path,
        default=ROOT / "vault" / "_shipping_finalists" / "manifest.json",
    )
    parser.add_argument(
        "--targets",
        type=Path,
        default=ROOT / "datasets" / "role_vocab" / "shipping_role_targets_v1.json",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "vault" / "_scorecards" / "role_vocab_scorecard.json",
    )
    args = parser.parse_args()

    if not args.manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {args.manifest}")
    if not args.targets.exists():
        raise FileNotFoundError(f"Targets not found: {args.targets}")

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    targets = json.loads(args.targets.read_text(encoding="utf-8")).get("targets", {})

    out = {
        "manifest": str(args.manifest),
        "targets": str(args.targets),
        "summary": [],
        "results": {},
    }

    for body_key, public_name in BODY_KEYS.items():
        target = targets.get(public_name)
        if target is None:
            raise KeyError(f"Missing role target for {public_name}")
        target_sig = target["signature"]
        rows = []

        for entry in manifest.get(body_key, []):
            candidate = Body.from_json(entry["path"])
            sig = body_signature(candidate)
            dist = role_distance(sig, target_sig)
            gate = shipping_gate(candidate, public_name)
            rows.append({
                "candidate_name": candidate.name,
                "candidate_path": entry["path"],
                "role_distance": dist,
                "shipping_gate_passed": gate["passed"],
                "shipping_gate_failures": gate["failures"],
            })

        rows.sort(key=lambda r: (
            0 if r["shipping_gate_passed"] else 1,
            r["role_distance"]["composite_distance"],
        ))
        chosen = rows[0] if rows else None
        out["results"][body_key] = {
            "public_name": public_name,
            "target_source": target["source_p2k"],
            "target_name": target["source_name"],
            "chosen": chosen,
            "candidates": rows,
        }
        if chosen is not None:
            out["summary"].append({
                "body": public_name,
                "candidate": chosen["candidate_name"],
                "shipping_gate_passed": chosen["shipping_gate_passed"],
                "role_composite_distance": chosen["role_distance"]["composite_distance"],
                "role_mismatch_rate": chosen["role_distance"]["role_mismatch_rate"],
                "role_mean_l2": chosen["role_distance"]["mean_l2"],
            })

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2), encoding="utf-8")
    print(args.out)


if __name__ == "__main__":
    main()
