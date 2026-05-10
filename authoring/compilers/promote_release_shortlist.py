#!/usr/bin/env python3
"""Promote a shortlist directory into a release-candidate package folder."""

from __future__ import annotations

import argparse
import json
import shutil
from pathlib import Path


def _load_shortlist(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _copy_entries(source_dir: Path, dest_dir: Path, entries: list[dict]) -> list[str]:
    dest_dir.mkdir(parents=True, exist_ok=True)
    copied: list[str] = []
    for e in entries:
        name = str(e["file"])
        src = source_dir / name
        dst = dest_dir / name
        if not src.exists():
            raise FileNotFoundError(f"missing shortlist file: {src}")
        shutil.copyfile(src, dst)
        copied.append(name)
    return copied


def main() -> int:
    ap = argparse.ArgumentParser(description="Promote shortlist to release candidates")
    ap.add_argument(
        "--shortlist-manifest",
        type=Path,
        required=True,
        help="path to manifest.shortlist.json",
    )
    ap.add_argument(
        "--shortlist-dir",
        type=Path,
        required=True,
        help="directory containing shortlist cartridge files",
    )
    ap.add_argument(
        "--release-dir",
        type=Path,
        required=True,
        help="destination release-candidate folder",
    )
    ap.add_argument(
        "--release-name",
        type=str,
        default="emu_style_release_candidate",
        help="identifier recorded in release manifest",
    )
    args = ap.parse_args()

    shortlist = _load_shortlist(args.shortlist_manifest)
    entries = list(shortlist.get("entries", []))
    if not entries:
        raise SystemExit("shortlist manifest has no entries")

    copied = _copy_entries(args.shortlist_dir, args.release_dir, entries)
    release_manifest = {
        "release_name": args.release_name,
        "source_shortlist_manifest": str(args.shortlist_manifest),
        "count": len(copied),
        "files": sorted(copied),
        "entries": entries,
    }
    (args.release_dir / "release_manifest.json").write_text(
        json.dumps(release_manifest, indent=2) + "\n",
        encoding="utf-8",
    )

    setlist = args.shortlist_dir / "SETLIST.md"
    if setlist.exists():
        shutil.copyfile(setlist, args.release_dir / "SETLIST.md")

    print(f"promoted {len(copied)} files to {args.release_dir}")
    print(args.release_dir / "release_manifest.json")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

