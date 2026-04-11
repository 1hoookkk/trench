"""Build role vocabulary from P2k corpus (P2k_000..P2k_032)."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.role_vocab import body_signature, summarize_role_vocabulary


def main() -> None:
    parser = argparse.ArgumentParser(description="Build P2k role vocabulary.")
    parser.add_argument(
        "--p2k-dir",
        type=Path,
        default=ROOT / "cartridges" / "p2k",
    )
    parser.add_argument(
        "--names",
        type=Path,
        default=ROOT / "datasets" / "p2k_filter_names.json",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=ROOT / "datasets" / "role_vocab" / "p2k_role_vocabulary_v1.json",
    )
    args = parser.parse_args()

    if not args.p2k_dir.exists():
        raise FileNotFoundError(f"P2k dir not found: {args.p2k_dir}")
    if not args.names.exists():
        raise FileNotFoundError(f"Names file not found: {args.names}")

    names_obj = json.loads(args.names.read_text(encoding="utf-8"))
    name_map = {int(k): v for k, v in names_obj.get("names", {}).items()}

    rows: list[tuple[str, dict]] = []
    body_rows = []
    for idx in range(33):
        path = args.p2k_dir / f"P2k_{idx:03d}.json"
        if not path.exists():
            continue
        body = Body.from_json(str(path))
        signature = body_signature(body)
        display_name = name_map.get(idx, body.name)
        key = f"P2k_{idx:03d}"
        rows.append((key, signature))
        body_rows.append({
            "index": idx,
            "key": key,
            "name": display_name,
            "path": str(path),
            "signature": signature,
        })

    vocab = summarize_role_vocabulary(rows)
    payload = {
        "version": "p2k-role-vocab-v1",
        "source": str(args.p2k_dir),
        "body_count": len(body_rows),
        "bodies": body_rows,
        "role_prototypes": vocab["role_prototypes"],
        "stage_role_distribution": vocab["stage_role_distribution"],
        "body_role_sequences": vocab["body_role_sequences"],
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(args.out)
    print(f"body_count={len(body_rows)}")
    print(f"roles={len(vocab['role_prototypes'])}")


if __name__ == "__main__":
    main()
