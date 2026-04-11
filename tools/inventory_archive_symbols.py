"""Inventory archive symbols against the live Trench workspace.

Outputs:
- JSON machine-readable manifest
- Markdown summary with key symbol status table
"""
from __future__ import annotations

import argparse
import json
import re
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_ARCHIVE_ROOT = Path("C:/Users/hooki/trenchwork_clean/archive")


RUST_PUB_RE = re.compile(
    r"^\s*pub\s+(?:\([^)]*\)\s+)?(?:(?:async|const|unsafe|extern)\s+)*"
    r"(?:fn|struct|enum|trait|type|mod|static)\s+([A-Za-z_][A-Za-z0-9_]*)"
)
RUST_FN_RE = re.compile(r"^\s*(?:async\s+)?fn\s+([A-Za-z_][A-Za-z0-9_]*)")
RUST_TYPE_RE = re.compile(r"^\s*(?:struct|enum|trait|type)\s+([A-Za-z_][A-Za-z0-9_]*)")
PY_DEF_RE = re.compile(r"^\s*def\s+([A-Za-z_][A-Za-z0-9_]*)\s*\(")
PY_CLASS_RE = re.compile(r"^\s*class\s+([A-Za-z_][A-Za-z0-9_]*)\b")

HERITAGE_PATTERNS = (
    "p2k",
    "morpheus",
    "emu",
    "historical",
    "heritage",
    "rom",
    "morphdesigner",
)
HERITAGE_TOKEN_RE = re.compile(
    r"(^|[^a-z0-9])(" + "|".join(HERITAGE_PATTERNS) + r")([^a-z0-9]|$)"
)
ADAPTATION_MODULE_HINTS = (
    "trench-atlas/src",
    "trench-station/src",
)

FOCUS_SYMBOLS = (
    "InterferenceLink",
    "ClusterPlan",
    "Role",
    "owned_gaps",
    "seed_relative_zero",
    "pressure_law",
    "prominence",
    "topology_anchor_theft",
    "topology_interference_drift",
    "topology_cluster_leakage",
    "topology_bypass_wake",
    "Material",
    "Motor",
    "Volatility",
    "Curation",
    "House",
)


@dataclass(frozen=True)
class SymbolRow:
    name: str
    kind: str
    visibility: str
    source_relpath: str
    line: int
    classification: str
    alive_in_shipping: bool
    dead_heritage: bool
    portable_to_current: bool
    needs_adaptation: bool
    live_hits: list[str]


def _iter_files(root: Path, suffixes: tuple[str, ...]) -> list[Path]:
    if not root.exists():
        return []
    files: list[Path] = []
    for suffix in suffixes:
        files.extend(root.rglob(f"*{suffix}"))
    return files


def _parse_archive_rust_symbols(base_root: Path, file_path: Path) -> list[dict]:
    symbols: list[dict] = []
    try:
        text = file_path.read_text(encoding="utf-8")
    except UnicodeDecodeError:
        text = file_path.read_text(encoding="latin-1")

    for idx, line in enumerate(text.splitlines(), start=1):
        pub_match = RUST_PUB_RE.match(line)
        if pub_match:
            name = pub_match.group(1)
            kind = "pub_item"
            if " fn " in f" {line} " or line.strip().startswith("pub fn"):
                kind = "pub_fn"
            elif " struct " in f" {line} ":
                kind = "pub_struct"
            elif " enum " in f" {line} ":
                kind = "pub_enum"
            symbols.append(
                {
                    "name": name,
                    "kind": kind,
                    "visibility": "pub",
                    "source_relpath": file_path.relative_to(base_root).as_posix(),
                    "line": idx,
                }
            )
            continue

        fn_match = RUST_FN_RE.match(line)
        if fn_match:
            symbols.append(
                {
                    "name": fn_match.group(1),
                    "kind": "fn",
                    "visibility": "private",
                    "source_relpath": file_path.relative_to(base_root).as_posix(),
                    "line": idx,
                }
            )
            continue

        type_match = RUST_TYPE_RE.match(line)
        if type_match:
            symbols.append(
                {
                    "name": type_match.group(1),
                    "kind": "type_item",
                    "visibility": "private",
                    "source_relpath": file_path.relative_to(base_root).as_posix(),
                    "line": idx,
                }
            )
    return symbols


def _collect_live_symbols(live_root: Path) -> tuple[set[str], dict[str, list[str]]]:
    live_dirs = [
        live_root / "trench-core",
        live_root / "trench-plugin",
        live_root / "trench-codec",
        live_root / "xtask",
        live_root / "pyruntime",
        live_root / "tools",
    ]
    symbols: set[str] = set()
    locations: dict[str, list[str]] = {}

    rust_files: list[Path] = []
    py_files: list[Path] = []
    for d in live_dirs:
        rust_files.extend(_iter_files(d, (".rs",)))
        py_files.extend(_iter_files(d, (".py",)))

    for file_path in rust_files:
        try:
            text = file_path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = file_path.read_text(encoding="latin-1")
        rel = file_path.relative_to(live_root).as_posix()
        for idx, line in enumerate(text.splitlines(), start=1):
            for rx in (RUST_PUB_RE, RUST_FN_RE, RUST_TYPE_RE):
                m = rx.match(line)
                if not m:
                    continue
                name = m.group(1)
                symbols.add(name)
                key = f"{rel}:{idx}"
                locations.setdefault(name, [])
                if len(locations[name]) < 8:
                    locations[name].append(key)
                break

    for file_path in py_files:
        try:
            text = file_path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            text = file_path.read_text(encoding="latin-1")
        rel = file_path.relative_to(live_root).as_posix()
        for idx, line in enumerate(text.splitlines(), start=1):
            for rx in (PY_DEF_RE, PY_CLASS_RE):
                m = rx.match(line)
                if not m:
                    continue
                name = m.group(1)
                symbols.add(name)
                key = f"{rel}:{idx}"
                locations.setdefault(name, [])
                if len(locations[name]) < 8:
                    locations[name].append(key)
                break
    return symbols, locations


def _classify_symbol(
    name: str,
    source_relpath: str,
    live_symbols: set[str],
) -> tuple[str, bool, bool, bool, bool]:
    lower_name = name.lower()
    lower_path = source_relpath.lower()

    alive = name in live_symbols
    dead_heritage = bool(HERITAGE_TOKEN_RE.search(lower_name)) or bool(
        HERITAGE_TOKEN_RE.search(lower_path)
    )
    needs_adaptation = any(h in lower_path for h in ADAPTATION_MODULE_HINTS)
    portable = (not alive) and (not dead_heritage) and (not needs_adaptation)

    if alive:
        classification = "alive_in_shipping"
    elif dead_heritage:
        classification = "dead_heritage"
    elif needs_adaptation:
        classification = "needs_adaptation"
    else:
        classification = "portable_to_current"

    return classification, alive, dead_heritage, portable, needs_adaptation


def _build_markdown(
    archive_root: Path,
    live_root: Path,
    rows: list[SymbolRow],
    focus_payload: list[dict],
    summary: dict[str, int],
) -> str:
    lines: list[str] = []
    lines.append("# Archive Symbol Inventory")
    lines.append("")
    lines.append(f"- Archive root: `{archive_root.as_posix()}`")
    lines.append(f"- Live root: `{live_root.as_posix()}`")
    lines.append(f"- Total symbols: `{summary['total_symbols']}`")
    lines.append(f"- `alive_in_shipping`: `{summary['alive_in_shipping']}`")
    lines.append(f"- `dead_heritage`: `{summary['dead_heritage']}`")
    lines.append(f"- `portable_to_current`: `{summary['portable_to_current']}`")
    lines.append(f"- `needs_adaptation`: `{summary['needs_adaptation']}`")
    lines.append("")
    lines.append("## Focus Symbol Status")
    lines.append("")
    lines.append("| Symbol | Status | Archive Source | Live Hits |")
    lines.append("|---|---|---|---|")
    for item in focus_payload:
        live_hits = ", ".join(item["live_hits"]) if item["live_hits"] else "-"
        lines.append(
            f"| {item['symbol']} | {item['status']} | "
            f"{item['archive_source'] or '-'} | {live_hits} |"
        )
    lines.append("")
    lines.append("## Portable Candidates (Top 50)")
    lines.append("")
    lines.append("| Symbol | Kind | Source |")
    lines.append("|---|---|---|")
    portable = [r for r in rows if r.classification == "portable_to_current"]
    for row in portable[:50]:
        lines.append(f"| {row.name} | {row.kind} | {row.source_relpath}:{row.line} |")
    lines.append("")
    lines.append("## Needs Adaptation (Top 50)")
    lines.append("")
    lines.append("| Symbol | Kind | Source |")
    lines.append("|---|---|---|")
    needs_adapt = [r for r in rows if r.classification == "needs_adaptation"]
    for row in needs_adapt[:50]:
        lines.append(f"| {row.name} | {row.kind} | {row.source_relpath}:{row.line} |")
    lines.append("")
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Cross-reference archive symbols against live shipping code."
    )
    parser.add_argument(
        "--archive-root",
        type=Path,
        default=DEFAULT_ARCHIVE_ROOT,
        help="Archive root containing pre-runtime/trench-atlas/trench-forge-legacy.",
    )
    parser.add_argument(
        "--live-root",
        type=Path,
        default=ROOT,
        help="Live repo root to compare against.",
    )
    parser.add_argument(
        "--out-json",
        type=Path,
        default=ROOT / "docs" / "archive_symbol_inventory_2026-04-09.json",
    )
    parser.add_argument(
        "--out-md",
        type=Path,
        default=ROOT / "docs" / "archive_symbol_inventory_2026-04-09.md",
    )
    args = parser.parse_args()

    archive_dirs = [
        args.archive_root / "pre-runtime" / "trench-forge" / "src",
        args.archive_root / "pre-runtime" / "trench-foundry" / "src",
        args.archive_root / "pre-runtime" / "trench-station" / "src",
        args.archive_root / "trench-atlas" / "src",
        args.archive_root / "trench-forge-legacy" / "src",
    ]
    existing_archive_dirs = [d for d in archive_dirs if d.exists()]
    if not existing_archive_dirs:
        raise FileNotFoundError(
            f"No archive source dirs found under {args.archive_root.as_posix()}"
        )

    live_symbols, live_locations = _collect_live_symbols(args.live_root)

    raw_symbols: list[dict] = []
    for d in existing_archive_dirs:
        for file_path in _iter_files(d, (".rs",)):
            raw_symbols.extend(_parse_archive_rust_symbols(args.archive_root, file_path))

    rows: list[SymbolRow] = []
    for sym in raw_symbols:
        (
            classification,
            alive,
            dead_heritage,
            portable,
            needs_adaptation,
        ) = _classify_symbol(sym["name"], sym["source_relpath"], live_symbols)
        rows.append(
            SymbolRow(
                name=sym["name"],
                kind=sym["kind"],
                visibility=sym["visibility"],
                source_relpath=sym["source_relpath"],
                line=sym["line"],
                classification=classification,
                alive_in_shipping=alive,
                dead_heritage=dead_heritage,
                portable_to_current=portable,
                needs_adaptation=needs_adaptation,
                live_hits=live_locations.get(sym["name"], []),
            )
        )

    rows_sorted = sorted(rows, key=lambda r: (r.classification, r.source_relpath, r.line))
    summary = {
        "total_symbols": len(rows_sorted),
        "alive_in_shipping": sum(r.alive_in_shipping for r in rows_sorted),
        "dead_heritage": sum(r.dead_heritage for r in rows_sorted),
        "portable_to_current": sum(r.portable_to_current for r in rows_sorted),
        "needs_adaptation": sum(r.needs_adaptation for r in rows_sorted),
    }

    focus_payload: list[dict] = []
    by_name: dict[str, list[SymbolRow]] = {}
    for row in rows_sorted:
        by_name.setdefault(row.name, []).append(row)

    for symbol in FOCUS_SYMBOLS:
        candidates = by_name.get(symbol, [])
        if candidates:
            row = candidates[0]
            focus_payload.append(
                {
                    "symbol": symbol,
                    "status": row.classification,
                    "archive_source": f"{row.source_relpath}:{row.line}",
                    "live_hits": row.live_hits,
                }
            )
            continue

        focus_payload.append(
            {
                "symbol": symbol,
                "status": "missing_in_archive_scan",
                "archive_source": None,
                "live_hits": live_locations.get(symbol, []),
            }
        )

    payload = {
        "generated_at": "2026-04-09",
        "archive_root": args.archive_root.as_posix(),
        "live_root": args.live_root.as_posix(),
        "summary": summary,
        "focus_symbols": focus_payload,
        "symbols": [
            {
                "name": r.name,
                "kind": r.kind,
                "visibility": r.visibility,
                "source_relpath": r.source_relpath,
                "line": r.line,
                "classification": r.classification,
                "alive_in_shipping": r.alive_in_shipping,
                "dead_heritage": r.dead_heritage,
                "portable_to_current": r.portable_to_current,
                "needs_adaptation": r.needs_adaptation,
                "live_hits": r.live_hits,
            }
            for r in rows_sorted
        ],
    }

    args.out_json.parent.mkdir(parents=True, exist_ok=True)
    args.out_md.parent.mkdir(parents=True, exist_ok=True)
    args.out_json.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    args.out_md.write_text(
        _build_markdown(args.archive_root, args.live_root, rows_sorted, focus_payload, summary),
        encoding="utf-8",
    )

    print(args.out_json.as_posix())
    print(args.out_md.as_posix())
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
