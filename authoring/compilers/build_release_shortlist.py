#!/usr/bin/env python3
"""Build a release shortlist from an emu_style_pack manifest.

Produces a smaller, balanced set for packaging/demo without manual file picking.
"""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


def _load_manifest(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _intensity_distance(value: float, target: float) -> float:
    return abs(value - target)


def _pick_entries(entries: list[dict], target_intensity: float, per_skin: int) -> list[dict]:
    # Prefer target intensity, then stronger presets, then category spread.
    buckets: dict[int, list[dict]] = {}
    for e in entries:
        buckets.setdefault(int(e["skin"]), []).append(e)

    picked: list[dict] = []
    for skin in sorted(buckets):
        skin_entries = buckets[skin]
        skin_entries.sort(
            key=lambda e: (
                _intensity_distance(float(e["intensity"]), target_intensity),
                -float(e.get("selection_score", 0.0)),
                str(e["category"]),
                str(e["target"]),
                str(e["recipe"]),
                str(e.get("variant", "")),
                str(e["file"]),
            )
        )
        # Keep one per recipe/variant family first to maximize variety.
        seen_signature: set[tuple[str, str]] = set()
        first_pass: list[dict] = []
        second_pass: list[dict] = []
        for e in skin_entries:
            signature = (str(e["recipe"]), str(e.get("variant", "")))
            if signature in seen_signature:
                second_pass.append(e)
                continue
            seen_signature.add(signature)
            first_pass.append(e)
        ordered = first_pass + second_pass
        picked.extend(ordered[:per_skin])
    return picked


def _write_setlist(path: Path, selected: list[dict], title: str) -> None:
    lines = [f"# {title}", ""]
    lines.append(f"- count: {len(selected)}")
    lines.append("")
    lines.append("| file | recipe | variant | skin | intensity | score | category | target |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for e in sorted(selected, key=lambda x: str(x["file"])):
        lines.append(
            f"| `{e['file']}` | `{e['recipe']}` | `{e.get('variant', 'core')}` | `{e['skin']}` | "
            f"`{e['intensity']}` | `{e.get('selection_score', 0.0)}` | `{e['category']}` | `{e['target']}` |"
        )
    lines.append("")
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    ap = argparse.ArgumentParser(description="Build shortlist from emu style manifest")
    ap.add_argument(
        "--manifest",
        type=Path,
        required=True,
        help="path to manifest.json from generate_emu_style_bank.py",
    )
    ap.add_argument(
        "--source-dir",
        type=Path,
        required=True,
        help="directory containing the source cartridge JSON files",
    )
    ap.add_argument(
        "--out-dir",
        type=Path,
        required=True,
        help="destination directory for shortlisted files",
    )
    ap.add_argument(
        "--target-intensity",
        type=float,
        default=1.0,
        help="preferred intensity to bias selection toward",
    )
    ap.add_argument(
        "--per-skin",
        type=int,
        default=6,
        help="number of files selected per skin",
    )
    ap.add_argument(
        "--title",
        type=str,
        default="E-mu Style Release Shortlist",
        help="title written into SETLIST.md",
    )
    args = ap.parse_args()

    manifest = _load_manifest(args.manifest)
    entries = list(manifest.get("entries", []))
    if not entries:
        raise SystemExit("manifest has no entries")
    if args.per_skin <= 0:
        raise SystemExit("per-skin must be > 0")

    selected = _pick_entries(entries, args.target_intensity, args.per_skin)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    for e in selected:
        src = args.source_dir / str(e["file"])
        dst = args.out_dir / str(e["file"])
        if not src.exists():
            raise SystemExit(f"missing source file: {src}")
        shutil.copyfile(src, dst)
        copied += 1

    shortlist_manifest = {
        "title": args.title,
        "count": copied,
        "target_intensity": args.target_intensity,
        "per_skin": args.per_skin,
        "entries": sorted(selected, key=lambda x: str(x["file"])),
    }
    (args.out_dir / "manifest.shortlist.json").write_text(
        json.dumps(shortlist_manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    _write_setlist(args.out_dir / "SETLIST.md", selected, args.title)

    print(f"shortlisted {copied} files into {args.out_dir}")
    print(args.out_dir / "SETLIST.md")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
