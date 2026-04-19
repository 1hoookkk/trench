"""Rapid audition logger — hot-swap queued pills into the plugin and log keystrokes.

Flow:
  1. Reads compiled-v1 pills from _queue/ in filename order.
  2. Writes each pill to TRENCH_LIVE_PATH (env var) or ~/Trench/trench_live.json.
     The running plugin polls that path and reloads on file change (~100 ms).
  3. Waits for a single keypress:
       1 = reject    2 = maybe    3 = keep    q/Esc = quit
  4. Appends one CSV row per verdict.

Usage (from repo root):
    python tools/workbench/verdict.py
    python tools/workbench/verdict.py --queue-dir tools/workbench/_queue --log verdicts.csv
    python tools/workbench/verdict.py --resume   # skip pills already in the log
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

QUEUE_DIR = Path(__file__).parent / "_queue"
CSV_PATH = Path(__file__).parent / "verdicts.csv"
CSV_FIELDS = ["timestamp_utc", "pill_id", "pill_name", "scope", "verdict", "file"]


def resolve_live_path() -> Path:
    env = os.environ.get("TRENCH_LIVE_PATH")
    if env:
        return Path(env)
    return Path.home() / "Trench" / "trench_live.json"


def get_keypress() -> str:
    if sys.platform == "win32":
        import msvcrt
        return msvcrt.getwch()   # single char, no ENTER needed
    else:
        import tty
        import termios
        fd = sys.stdin.fileno()
        old = termios.tcgetattr(fd)
        try:
            tty.setraw(fd)
            return sys.stdin.read(1)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old)


def pill_meta(pill: dict) -> tuple[str, str, str]:
    prov = pill.get("provenance", {})
    pill_id = prov.get("id") or pill.get("name", "?")
    name = prov.get("name") or pill.get("name", "?")
    scope = prov.get("scope") or pill.get("scope") or "?"
    return pill_id, name, scope


VERDICT_MAP = {"1": "reject", "2": "maybe", "3": "keep"}
QUIT_KEYS = {"q", "\x1b", "\x03"}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Audition queued pills in the live plugin and log verdicts.",
    )
    parser.add_argument("--queue-dir", type=Path, default=QUEUE_DIR,
                        help="Directory of compiled-v1 pills (default: _queue/).")
    parser.add_argument("--log", type=Path, default=CSV_PATH,
                        help="CSV verdict log (default: verdicts.csv).")
    parser.add_argument("--resume", action="store_true",
                        help="Skip pills already recorded in the log.")
    args = parser.parse_args(argv)

    live_path = resolve_live_path()
    queue_files = sorted(args.queue_dir.glob("*.json"))
    if not queue_files:
        print(f"no pills in {args.queue_dir} — run sampler.py first", file=sys.stderr)
        return 1

    already_logged: set[str] = set()
    if args.resume and args.log.exists():
        with open(args.log, newline="", encoding="utf-8") as f:
            for row in csv.DictReader(f):
                already_logged.add(row.get("file", ""))

    log_exists = args.log.exists()
    log_fh = open(args.log, "a", newline="", encoding="utf-8")
    writer = csv.DictWriter(log_fh, fieldnames=CSV_FIELDS)
    if not log_exists:
        writer.writeheader()

    print(f"live path : {live_path}")
    print(f"log       : {args.log}")
    print(f"queue     : {len(queue_files)} pills")
    print("controls  : 1=reject  2=maybe  3=keep  q=quit")
    print()

    kept = maybe = rejected = skipped = 0

    for i, fpath in enumerate(queue_files, 1):
        rel = str(fpath)
        if rel in already_logged:
            skipped += 1
            continue

        pill = json.loads(fpath.read_text(encoding="utf-8"))
        pill_id, name, scope = pill_meta(pill)

        shutil.copy2(fpath, live_path)

        print(f"[{i:3d}/{len(queue_files)}] {name:<40s} (scope={scope})  ", end="", flush=True)

        ch = get_keypress()

        if ch in QUIT_KEYS:
            print("(quit)")
            break

        verdict_str = VERDICT_MAP.get(ch)
        if verdict_str is None:
            print(f"(unknown key {ch!r} — skipping)")
            continue

        if verdict_str == "reject":
            rejected += 1
        elif verdict_str == "maybe":
            maybe += 1
        else:
            kept += 1

        print(verdict_str)

        writer.writerow({
            "timestamp_utc": datetime.now(timezone.utc).isoformat(),
            "pill_id": pill_id,
            "pill_name": name,
            "scope": scope,
            "verdict": verdict_str,
            "file": rel,
        })
        log_fh.flush()

    log_fh.close()
    print(f"\nsession: keep={kept}  maybe={maybe}  reject={rejected}  skipped={skipped}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
