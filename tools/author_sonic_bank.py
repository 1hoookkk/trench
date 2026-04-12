"""Authored shape bank: one body per sonic table target.

No RNG. No filtering. No bucketing. Walk the sonic tables and build a body
for every known target using the existing target.py builders + macro_compile.

Output: vault/_shapes/<category>/<key>.json
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.macro_compile import compile_body
from pyruntime.target import (
    build_vowel_target,
    build_nasal_target,
    build_landmark_target,
    build_bell_target,
    build_electronic_target,
    build_consonant_target,
    build_instrument_target,
    build_moog_target,
    build_singer_target,
    build_morph_target,
)

TABLES_PATH = ROOT / "docs" / "sonic_tables" / "tables.json"


def _slug(text: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in text).strip("_")


def _author_one(builder, key: str, name_prefix: str, out_dir: Path) -> tuple[bool, str]:
    try:
        spec = builder(f"{name_prefix}_{_slug(key)}", key)
        corners = compile_body(spec)
        body = Body(name=spec.name, corners=corners, boost=spec.boost)
        out_path = out_dir / f"{_slug(key)}.json"
        out_path.write_text(body.to_compiled_json(provenance="author-sonic-bank"), encoding="utf-8")
        return True, str(out_path)
    except Exception as e:
        return False, f"{key}: {e}"


def main():
    parser = argparse.ArgumentParser(description="Author one body per sonic table target.")
    parser.add_argument("--out-root", type=Path, default=ROOT / "vault" / "_shapes")
    parser.add_argument("--tables", type=Path, default=TABLES_PATH)
    parser.add_argument("--pairs", action="store_true",
                        help="Also author morph-pair bodies (vowel→vowel, vowel→nasal, etc).")
    args = parser.parse_args()

    tables = json.loads(args.tables.read_text(encoding="utf-8"))

    # Category → (keys, builder fn)
    categories: list[tuple[str, list[str], callable]] = []

    if "vowels" in tables:
        categories.append(("vowels", list(tables["vowels"].keys()), build_vowel_target))
    if "nasals" in tables:
        categories.append(("nasals", list(tables["nasals"].keys()), build_nasal_target))
    if "landmarks" in tables:
        names = [e["name"] for e in tables["landmarks"].get("entries", [])]
        categories.append(("landmarks", names, build_landmark_target))
    if "consonants" in tables:
        categories.append(("consonants", list(tables["consonants"].keys()), build_consonant_target))
    if "instruments" in tables:
        categories.append(("instruments", list(tables["instruments"].keys()), build_instrument_target))
    if "moog" in tables:
        categories.append(("moog", list(tables["moog"].keys()), build_moog_target))
    if "singer" in tables:
        categories.append(("singer", list(tables["singer"].keys()), build_singer_target))
    if "electronic" in tables:
        categories.append(("electronic", list(tables["electronic"].keys()), build_electronic_target))
    if "measured_bells" in tables:
        bells = tables["measured_bells"]
        bell_keys = [b.get("name") for b in bells if isinstance(b, dict) and b.get("name")]
        categories.append(("measured_bells", bell_keys, build_bell_target))

    total_ok = 0
    total_fail = 0
    fails: list[str] = []

    for cat_name, keys, builder in categories:
        cat_dir = args.out_root / cat_name
        cat_dir.mkdir(parents=True, exist_ok=True)
        for key in keys:
            ok, msg = _author_one(builder, key, cat_name, cat_dir)
            if ok:
                total_ok += 1
                print(f"  OK  {cat_name}/{key}")
            else:
                total_fail += 1
                fails.append(msg)
                print(f"  ERR {cat_name}/{msg}", file=sys.stderr)

    # Morph pairs: vowel↔vowel, vowel↔nasal, vowel↔consonant
    if args.pairs:
        pair_defs = [
            ("vowel_pairs",  "vowel", "vowel",     tables.get("vowels", {}),  tables.get("vowels", {})),
            ("vowel_nasal",  "vowel", "nasal",     tables.get("vowels", {}),  tables.get("nasals", {})),
            ("vowel_cons",   "vowel", "consonant", tables.get("vowels", {}),  tables.get("consonants", {})),
            ("vowel_inst",   "vowel", "instrument",tables.get("vowels", {}),  tables.get("instruments", {})),
            ("vowel_land",   "vowel", "landmark",  tables.get("vowels", {}),  {e["name"]: e for e in tables.get("landmarks", {}).get("entries", [])}),
        ]
        for cat_name, src_a, src_b, keys_a, keys_b in pair_defs:
            pair_dir = args.out_root / cat_name
            pair_dir.mkdir(parents=True, exist_ok=True)
            for ka in keys_a:
                for kb in keys_b:
                    if ka == kb:
                        continue
                    try:
                        name = f"{cat_name}_{_slug(ka)}_to_{_slug(kb)}"
                        spec = build_morph_target(name, src_a, ka, src_b, kb)
                        corners = compile_body(spec)
                        body = Body(name=name, corners=corners, boost=spec.boost)
                        fp = pair_dir / f"{_slug(ka)}_to_{_slug(kb)}.json"
                        fp.write_text(body.to_compiled_json(provenance="author-sonic-bank"), encoding="utf-8")
                        total_ok += 1
                        print(f"  OK  {cat_name}/{_slug(ka)}_to_{_slug(kb)}")
                    except Exception as e:
                        total_fail += 1
                        fails.append(f"{cat_name}/{ka}->{kb}: {e}")
                        print(f"  ERR {cat_name}/{ka}->{kb}: {e}", file=sys.stderr)

    # Write manifest
    manifest_path = args.out_root / "manifest.json"
    manifest_path.write_text(json.dumps({
        "authored": total_ok,
        "failed": total_fail,
        "fails": fails,
    }, indent=2), encoding="utf-8")

    print(f"\nAuthored: {total_ok}")
    print(f"Failed:   {total_fail}")
    print(f"Output:   {args.out_root}")


if __name__ == "__main__":
    main()
