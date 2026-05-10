"""Tiny HTTP server bridging the pole-authoring HTML to the standalone's
live-cartridge polling path.

Serves `index.html` on GET / and accepts POST /write with a compiled-v1
JSON body — writes that body to `~/trench-juce/trench_live.json` (the
path the JUCE plugin's `resolveLiveJsonPath` polls). Standalone hot-
reloads on file mtime change.

Stdlib only.

Usage:
    python tools/pole_authoring/server.py [--port 8765] [--target PATH]

Then open http://localhost:8765 in a browser.
"""
from __future__ import annotations

import argparse
import json
import os
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path

REPO = Path(__file__).resolve().parents[2]
DEFAULT_TARGET = Path(os.path.expanduser("~")) / "trench-juce" / "trench_live.json"
INDEX_HTML = Path(__file__).parent / "index.html"
CORPUS_DIR = REPO / "vault" / "_shapes"


_CORNER_SUFFIXES = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")


def _scan_templates() -> dict[str, dict[str, Path]]:
    """Walk vault/_shapes and group files by template stem (source + name).

    Returns {full_stem: {corner_label: path}}. A template is usable if it has
    at least M0_Q0 and M100_Q0; Q100 corners are optional (CORPUS only has
    Q0, so those get mirrored on synth; P2K has the full 4-corner set).
    """
    out: dict[str, dict[str, Path]] = {}
    if not CORPUS_DIR.is_dir():
        return out
    for p in CORPUS_DIR.glob("*.json"):
        if p.name == "shapes_index.json":
            continue
        stem = p.stem  # "CORPUS_AC_Crazy_Rez_M0_Q0"
        for suf in _CORNER_SUFFIXES:
            tag = f"_{suf}"
            if stem.endswith(tag):
                tstem = stem[: -len(tag)]
                out.setdefault(tstem, {})[suf] = p
                break
    return {s: c for s, c in out.items() if "M0_Q0" in c and "M100_Q0" in c}


def _corpus_pairs() -> list[dict]:
    """Return list of {name, has_full_q} sorted alphabetically (case-insensitive)."""
    scan = _scan_templates()
    items = []
    for stem, corners in scan.items():
        items.append({
            "name": stem,
            "has_full_q": "M0_Q100" in corners and "M100_Q100" in corners,
        })
    items.sort(key=lambda x: x["name"].lower())
    return items


def _synth_pair(stem: str) -> dict | None:
    corners = _scan_templates().get(stem)
    if corners is None:
        return None
    loaded: dict[str, list] = {}
    for suf, path in corners.items():
        loaded[suf] = json.loads(path.read_text(encoding="utf-8")).get("stages", [])
    # Mirror Q-row from M-row when genuine Q100 data is absent.
    m0_q100 = loaded.get("M0_Q100", loaded["M0_Q0"])
    m100_q100 = loaded.get("M100_Q100", loaded["M100_Q0"])
    return {
        "format": "compiled-v1",
        "name": stem,
        "provenance": f"vault/_shapes pair ({stem}; Q100={'native' if 'M0_Q100' in corners else 'mirrored'})",
        "sampleRate": 39062.5,
        "stages": 12,
        "keyframes": [
            {"label": "M0_Q0",     "morph": 0, "q": 0, "stages": loaded["M0_Q0"]},
            {"label": "M100_Q0",   "morph": 1, "q": 0, "stages": loaded["M100_Q0"]},
            {"label": "M0_Q100",   "morph": 0, "q": 1, "stages": m0_q100},
            {"label": "M100_Q100", "morph": 1, "q": 1, "stages": m100_q100},
        ],
    }


class Handler(BaseHTTPRequestHandler):
    target_path: Path = DEFAULT_TARGET

    def _send(self, code: int, body: bytes, ctype: str = "text/plain") -> None:
        self.send_response(code)
        self.send_header("Content-Type", ctype)
        self.send_header("Content-Length", str(len(body)))
        self.send_header("Access-Control-Allow-Origin", "*")
        self.end_headers()
        self.wfile.write(body)

    def do_GET(self) -> None:
        if self.path in ("/", "/index.html"):
            html = INDEX_HTML.read_bytes()
            self._send(200, html, "text/html; charset=utf-8")
            return
        if self.path == "/target":
            payload = json.dumps({"target": str(self.target_path)}).encode("utf-8")
            self._send(200, payload, "application/json")
            return
        if self.path == "/corpus/list":
            payload = json.dumps({"pairs": _corpus_pairs()}).encode("utf-8")
            self._send(200, payload, "application/json")
            return
        if self.path.startswith("/corpus/pair?"):
            from urllib.parse import urlparse, parse_qs
            q = parse_qs(urlparse(self.path).query)
            names = q.get("name", [])
            if not names:
                self._send(400, b"missing ?name=")
                return
            synth = _synth_pair(names[0])
            if synth is None:
                self._send(404, b"pair not found")
                return
            self._send(200, json.dumps(synth).encode("utf-8"), "application/json")
            return
        self._send(404, b"not found")

    def do_POST(self) -> None:
        if self.path != "/write":
            self._send(404, b"not found")
            return
        length = int(self.headers.get("Content-Length", "0"))
        raw = self.rfile.read(length).decode("utf-8")
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError as e:
            self._send(400, f"invalid JSON: {e}".encode("utf-8"))
            return
        try:
            self.target_path.parent.mkdir(parents=True, exist_ok=True)
            self.target_path.write_text(json.dumps(parsed, indent=2), encoding="utf-8")
        except OSError as e:
            self._send(500, f"write failed: {e}".encode("utf-8"))
            return
        msg = f"wrote {self.target_path} ({len(raw)} bytes)"
        self._send(200, msg.encode("utf-8"))

    def do_OPTIONS(self) -> None:
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.end_headers()

    def log_message(self, fmt: str, *args) -> None:  # noqa: D401
        # Compact log: METHOD path → status
        print(f"[pole-server] {self.address_string()} - {fmt % args}")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--port", type=int, default=8765)
    parser.add_argument("--target", type=Path, default=DEFAULT_TARGET,
                        help="trench_live.json path the standalone polls")
    args = parser.parse_args(argv)
    Handler.target_path = args.target

    server = ThreadingHTTPServer(("127.0.0.1", args.port), Handler)
    print(f"pole-authoring server: http://localhost:{args.port}")
    print(f"  writes → {args.target}")
    print("  ctrl-c to stop")
    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nstopped")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
