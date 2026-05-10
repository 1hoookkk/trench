#!/usr/bin/env python3
"""Forge write wrapper — compile raw-stage-v1 and atomically land it in the
trench-forge-dev watch slot.

Usage:
    python tools/forge_write.py -              # read raw-stage-v1 from stdin
    python tools/forge_write.py raw.json       # read from a file

The wrapper invokes `tools/compile_raw.py` with `--out <watch>.tmp`, then
`os.replace()`s the tmp file onto the watch path so the plugin's notify
watcher only ever sees a fully-written file. Same-directory rename keeps
the operation atomic on Windows (NTFS) and POSIX.

Watch slot: `$TRENCH_FORGE_WATCH`, default `forge-watch/active.json`.
The directory is created if missing.

Failure mode:
    - compile_raw.py exit != 0 → stderr is written to
      `<watch_dir>/last_error.log`, `.tmp` is unlinked, wrapper exits 1.
    - On success, `last_error.log` is overwritten with `ok`.
"""

from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

DEFAULT_WATCH = "forge-watch/active.json"


def watch_path() -> Path:
    return Path(os.environ.get("TRENCH_FORGE_WATCH", DEFAULT_WATCH))


def main(argv: list[str]) -> int:
    if len(argv) != 1:
        print("usage: forge_write.py (- | path/to/raw.json)", file=sys.stderr)
        return 2

    src = argv[0]
    repo_root = Path(__file__).resolve().parent.parent
    compile_raw = repo_root / "tools" / "compile_raw.py"

    out_path = watch_path()
    if not out_path.is_absolute():
        out_path = repo_root / out_path
    watch_dir = out_path.parent
    watch_dir.mkdir(parents=True, exist_ok=True)
    tmp_path = out_path.with_suffix(out_path.suffix + ".tmp")
    err_log = watch_dir / "last_error.log"

    stdin_bytes: bytes | None = None
    if src == "-":
        stdin_bytes = sys.stdin.buffer.read()

    proc = subprocess.run(
        [sys.executable, str(compile_raw), src, "--out", str(tmp_path)],
        input=stdin_bytes,
        capture_output=True,
    )

    if proc.returncode != 0:
        err = proc.stderr.decode("utf-8", errors="replace").rstrip()
        err_log.write_text(
            f"compile failed (exit {proc.returncode})\nsrc: {src}\n\n{err}\n"
        )
        try:
            tmp_path.unlink()
        except FileNotFoundError:
            pass
        print(err_log.read_text(), file=sys.stderr)
        return 1

    os.replace(tmp_path, out_path)
    err_log.write_text("ok\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
