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

# ---------------------------------------------------------------------------
# Pill labels — 1-4 character display names for a Looperator-style grid cell.
#
# Policy:
#   - Vowels use IPA-adjacent 2-char codes (AH, EE, OO, ...).
#   - Consonants and nasals are 1-2 chars (P, S, SH, M, N).
#   - Heritage phonemes use a 3-char H## shorthand (H04 for PH_004).
#   - Everything else gets a 2-4 char mnemonic picked per token.
#
# If a token isn't in this table the bake falls back to an uppercased
# version of its key truncated to 4 chars — never silently fails, but
# flags the token in stdout so we can decide whether to add a proper
# pill.
# ---------------------------------------------------------------------------
PILL_LABELS: dict[str, str] = {
    # consonants
    "consonants.p": "P",
    "consonants.s": "S",
    "consonants.sh": "SH",
    "consonants.t": "T",
    # electronic
    "electronic.acid": "ACID",
    "electronic.telephone": "PHN",
    # heritage phonemes
    "heritage.ph_004": "H04",
    "heritage.ph_009": "H09",
    "heritage.ph_011": "H11",
    "heritage.ph_013": "H13",
    "heritage.ph_021": "H21",
    "heritage.ph_023": "H23",
    "heritage.ph_028": "H28",
    # instruments
    "instruments.bassoon": "BSN",
    "instruments.clarinet": "CLR",
    "instruments.flute": "FLT",
    "instruments.french_horn": "FH",
    "instruments.oboe": "OB",
    "instruments.trumpet": "TRP",
    "instruments.tuba": "TBA",
    # landmarks
    "landmarks.air_band": "AIR",
    "landmarks.chest_resonance": "CHST",
    "landmarks.ear_canal_resonance": "ECR",
    "landmarks.presence_peak": "PRS",
    "landmarks.room_mode_medium": "RMM",
    "landmarks.room_mode_small": "RMS",
    "landmarks.sibilance_band": "SIB",
    "landmarks.telephone_band_high": "THI",
    "landmarks.telephone_band_low": "TLO",
    "landmarks.vocal_tract_f1_floor": "F1F",
    # measured bells
    "measured_bells.berlin_freedom": "BRL",
    "measured_bells.freiburg_hosanna": "FRB",
    "measured_bells.st_mary_le_tower": "STM",
    "measured_bells.stretched_treble": "STR",
    # moog
    "moog.moog_bass_compensated": "MGB",
    "moog.moog_self_oscillating": "MGO",
    "moog.moog_uncompensated": "MGU",
    # nasals
    "nasals.nasal_m": "M",
    "nasals.nasal_n": "N",
    # singer
    "singer.singers_formant": "SF",
    # vowels — IPA-style 2-char codes
    "vowels.ae": "AE",
    "vowels.ah": "AH",
    "vowels.aw": "AW",
    "vowels.ee": "EE",
    "vowels.eh": "EH",
    "vowels.er": "ER",
    "vowels.ih": "IH",
    "vowels.oo": "OO",
    "vowels.schwa": "SC",
    "vowels.uh": "UH",
}


def pill_label_for(token_id: str, key: str) -> tuple[str, bool]:
    """Return (pill_label, was_explicit). If the token isn't in the policy
    table, fall back to uppercased key truncated to 4 chars so nothing
    silently lands without a label — but flag it so the caller can add
    an explicit entry next bake."""
    if token_id in PILL_LABELS:
        return PILL_LABELS[token_id], True
    return key.upper()[:4], False


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
    unlabeled: list[str] = []
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

        pill, explicit = pill_label_for(token_id, key)
        if not explicit:
            unlabeled.append(f"{token_id} → fallback pill '{pill}'")

        manifest_entries.append(
            {
                "id": token_id,
                "pill": pill,
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
    if unlabeled:
        print(f"{len(unlabeled)} tokens fell back to auto pill labels:")
        for u in unlabeled:
            print(f"  - {u}")
    if missing:
        print(f"{len(missing)} inventory entries skipped (shape missing):")
        for m in missing:
            print(f"  - {m}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
