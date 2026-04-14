"""Bake phoneme pills from the unified token inventory → cartridges/engine/.

Reads `vault/_phonemes/token_inventory_unified_v2.json`, copies each token's
referenced shape file from `vault/_shapes/<category>/<key>.json` to
`cartridges/engine/<category>/<key>.json`, and writes a single manifest at
`cartridges/engine/manifest.json` mapping short labels to cartridge paths.

No re-authoring, no schema transformation. The shape files under
`vault/_shapes/` are already `compiled-v1` cartridges that load directly
through `trench_core::Cartridge::from_json`. This tool just lifts the ones
the canonical inventory names into the shipping directory, gives them a
stable path, and produces a Looperator-pill-friendly manifest.

Stdlib only. Runs in any clean venv. Re-runnable — each run is a clean
rebuild of `cartridges/engine/` from source.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
INVENTORY = REPO / "vault" / "_phonemes" / "token_inventory_unified_v2.json"
SHAPES = REPO / "vault" / "_shapes"
ENGINE = REPO / "cartridges" / "engine"


def main() -> int:
    if not INVENTORY.exists():
        raise SystemExit(f"missing inventory: {INVENTORY}")
    if not SHAPES.is_dir():
        raise SystemExit(f"missing shape bank: {SHAPES}")

    inventory = json.loads(INVENTORY.read_text(encoding="utf-8"))
    tokens = inventory.get("tokens", {})
    if not tokens:
        raise SystemExit(f"empty inventory: {INVENTORY}")

    # Clean rebuild — the engine dir is derived, so we own it fully.
    if ENGINE.exists():
        shutil.rmtree(ENGINE)
    ENGINE.mkdir(parents=True)

    manifest_entries: list[dict] = []
    missing: list[str] = []
    loaded = 0

    for token_id, meta in tokens.items():
        category = meta.get("category", "")
        key = meta.get("key", "")
        label = meta.get("label", token_id)
        if not category or not key:
            missing.append(f"{token_id} (no category/key)")
            continue

        src = SHAPES / category / f"{key}.json"
        if not src.exists():
            missing.append(f"{token_id} → {src.relative_to(REPO)}")
            continue

        shape = json.loads(src.read_text(encoding="utf-8"))
        dst_dir = ENGINE / category
        dst_dir.mkdir(parents=True, exist_ok=True)
        dst = dst_dir / f"{key}.json"
        dst.write_text(json.dumps(shape, indent=2) + "\n", encoding="utf-8")

        manifest_entries.append(
            {
                "id": token_id,
                "label": label,
                "category": category,
                "key": key,
                "path": str(dst.relative_to(REPO)).replace("\\", "/"),
            }
        )
        loaded += 1

    manifest = {
        "format": "phoneme-pill-manifest-v1",
        "generated_from": str(INVENTORY.relative_to(REPO)).replace("\\", "/"),
        "token_count": loaded,
        "tokens": manifest_entries,
    }
    (ENGINE / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    print(f"baked {loaded} phoneme pills → {ENGINE.relative_to(REPO)}/")
    print(f"wrote manifest → {(ENGINE / 'manifest.json').relative_to(REPO)}")
    if missing:
        print(f"{len(missing)} inventory entries skipped (shape missing):")
        for m in missing:
            print(f"  - {m}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
