from __future__ import annotations

import datetime
import json
import os
import pathlib
from typing import Any


_PRUNE_DIR_PREFIXES: tuple[str, ...] = (
    # Repo meta / caches
    ".git",
    ".agents",
    ".claude",
    ".codex",
    ".pytest_cache",
    ".vscode",
    "__pycache__",
    ".venv",
    "node_modules",
    # Big/non-canonical outputs
    "_GRAVEYARD",
    "target",
    "target-codex-vizia",
    "dev/tmp",
    "vault/_shapes",
    "vault/_profiles",
    "pyruntime/.terrain_cache",
    # JUCE build junk (requested)
    "trench-juce/plugin/build",
    "trench-juce/plugin/Builds",
    "trench-juce/plugin/JUCE",
    "trench-juce/plugin/.overnight_cargo_target",
    # MCP servers can include huge deps; not canonical truth.
    "mcp-servers-jos",
)


def _iso_now() -> str:
    # No timezone assumptions; keep it local-time ISO for easy reading.
    return datetime.datetime.now().isoformat(timespec="seconds")


def _read_json(path: pathlib.Path) -> dict[str, Any] | list[Any] | None:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError:
        return None


def _exists(rel: str) -> bool:
    return (REPO_ROOT / rel).exists()


def _posix_rel(path_abs: pathlib.Path) -> str:
    return str(path_abs.resolve().relative_to(REPO_ROOT)).replace("\\", "/")


def _should_prune_dir(rel_posix: str) -> bool:
    for prefix in _PRUNE_DIR_PREFIXES:
        if rel_posix == prefix or rel_posix.startswith(prefix + "/"):
            return True
    return False


def _is_stub(path: pathlib.Path) -> bool:
    """Heuristic: stubs must loudly declare themselves + point to canonical."""

    try:
        head = path.read_text(encoding="utf-8", errors="replace").splitlines()[:30]
    except Exception:
        return False

    joined = "\n".join(head)
    return ("GENERATED STUB" in joined) and ("Canonical:" in joined)


def _scan_truth_files() -> dict[str, Any]:
    """Find competing doctrine/law files and ensure duplicates are stubs.

    This is the mechanical enforcement that prevents “start directory decides doctrine.”
    The scan is pruned aggressively so it stays fast.
    """

    rules = {
        "SPEC.md": {"canonical": "SPEC.md", "allowed_non_stub": []},
        "SHIPPING.md": {"canonical": "SHIPPING.md", "allowed_non_stub": []},
        "AGENTS.md": {"canonical": "AGENTS.md", "allowed_non_stub": []},
        # CLAUDE.md is special: root is canonical doctrine, but subtrees may have scoped rules.
        "CLAUDE.md": {"canonical": "CLAUDE.md", "allowed_non_stub": ["pyruntime/CLAUDE.md", "trench-core/CLAUDE.md"]},
        # REPO_TRUTH is canonical at root; duplicates are not expected.
        "REPO_TRUTH.md": {"canonical": "REPO_TRUTH.md", "allowed_non_stub": []},
        "REPO_TRUTH.json": {"canonical": "REPO_TRUTH.json", "allowed_non_stub": []},
    }

    found: dict[str, list[dict[str, Any]]] = {k: [] for k in rules.keys()}
    violations: list[dict[str, Any]] = []

    for dirpath, dirnames, filenames in os.walk(REPO_ROOT):
        dir_abs = pathlib.Path(dirpath)
        try:
            rel_dir = dir_abs.resolve().relative_to(REPO_ROOT)
        except Exception:
            # Should never happen; fail closed (don’t walk unknown roots).
            dirnames[:] = []
            continue

        # Prune dirs in-place (keeps scan fast).
        for d in list(dirnames):
            rel_child = (rel_dir / d)
            rel_child_posix = str(rel_child).replace("\\", "/")
            if _should_prune_dir(rel_child_posix):
                dirnames.remove(d)
                continue

            # Do not descend into symlinked dirs/junctions; they explode the scan surface.
            try:
                if (dir_abs / d).is_symlink():
                    dirnames.remove(d)
            except OSError:
                # If we can’t stat it reliably, don’t descend.
                dirnames.remove(d)

        for fn in filenames:
            if fn not in rules:
                continue

            p = dir_abs / fn
            try:
                rel_path = _posix_rel(p)
                size = p.stat().st_size
            except OSError:
                rel_path = str(p)
                size = None

            canonical = rules[fn]["canonical"]
            allowed = set(rules[fn].get("allowed_non_stub") or [])

            entry: dict[str, Any] = {
                "path": rel_path,
                "size": size,
                "kind": "unknown",
                "stub_ok": None,
                "ok": True,
            }

            if rel_path == canonical:
                entry["kind"] = "canonical"
            elif rel_path in allowed:
                entry["kind"] = "allowed_non_stub"
            else:
                entry["kind"] = "stub"
                entry["stub_ok"] = _is_stub(p)
                entry["ok"] = bool(entry["stub_ok"])
                if not entry["ok"]:
                    violations.append(
                        {
                            "file": fn,
                            "path": rel_path,
                            "reason": "duplicate truth file is not a GENERATED STUB",
                        }
                    )

            found[fn].append(entry)

    # Sort for stability.
    for fn in found:
        found[fn].sort(key=lambda e: (e.get("kind") != "canonical", str(e.get("path"))))

    return {
        "generated_at": _iso_now(),
        "prune_dir_prefixes": list(_PRUNE_DIR_PREFIXES),
        "rules": rules,
        "found": found,
        "violations": violations,
        "ok": len(violations) == 0,
    }


def _shape_counts() -> dict[str, Any]:
    manifest_unified = REPO_ROOT / "vault/_shapes/manifest_unified.json"
    if (data := _read_json(manifest_unified)) and isinstance(data, dict):
        return {
            "source": "vault/_shapes/manifest_unified.json",
            "counts": data.get("counts", {}),
        }

    manifest = REPO_ROOT / "vault/_shapes/manifest.json"
    if (data := _read_json(manifest)) and isinstance(data, dict):
        authored = data.get("authored")
        failed = data.get("failed")
        return {
            "source": "vault/_shapes/manifest.json",
            "counts": {"total_bodies": authored, "failed": failed},
        }

    return {"source": None, "counts": {}}


def _p2k_counts() -> dict[str, Any]:
    inv_v2 = REPO_ROOT / "vault/_phonemes/p2k_phoneme_inventory_v2.json"
    aliases_v2 = REPO_ROOT / "vault/_phonemes/p2k_phoneme_aliases_v2.json"

    clusters = None
    unmapped = None
    heritage_fills = None

    if (data := _read_json(inv_v2)) and isinstance(data, dict):
        clusters = data.get("clusters")
        unmapped_list = data.get("unmapped_clusters")
        if isinstance(unmapped_list, list):
            unmapped = len(unmapped_list)

    if (data := _read_json(aliases_v2)) and isinstance(data, dict):
        aliases = data.get("aliases", {})
        if isinstance(aliases, dict):
            # Count clusters filled by heritage.* tokens.
            heritage_fills = 0
            for _, v in aliases.items():
                ref = None
                if isinstance(v, dict):
                    ref = v.get("alias_tables_ref")
                if isinstance(ref, str) and ref.startswith("heritage."):
                    heritage_fills += 1

            if clusters is None:
                clusters = len(aliases)
            if unmapped is None:
                unmapped = sum(1 for _, v in aliases.items() if not (isinstance(v, dict) and v.get("alias_tables_ref")))

    return {
        "clusters": clusters,
        "unmapped": unmapped,
        "heritage_fills": heritage_fills,
    }


def _external_workspaces() -> list[dict[str, Any]]:
    """Summarize other local roots that commonly get mixed with this repo.

    These are *not* canonical for this repo. We surface them so operators and agents don't
    accidentally treat them as part of the shipping codebase.
    """

    home = pathlib.Path.home()
    candidates: list[tuple[str, str]] = [
        ("trench_re_vault", str(home / "trench_re_vault")),
        ("trenchwork_clean", str(home / "trenchwork_clean")),
    ]

    out: list[dict[str, Any]] = []
    for name, path_abs in candidates:
        root = pathlib.Path(path_abs)
        exists = root.exists()
        entry: dict[str, Any] = {
            "name": name,
            "path_abs": path_abs,
            "exists": exists,
        }
        if not exists:
            out.append(entry)
            continue

        entry["git_repo"] = (root / ".git").exists()
        for f in ("AGENTS.md", "CLAUDE.md", "handoff.json"):
            entry[f"has_{f.lower().replace('.', '_')}"] = (root / f).exists()

        # Cheap "what is this" hint: first line of AGENTS.md if present.
        agents_path = root / "AGENTS.md"
        if agents_path.exists():
            try:
                first_line = agents_path.read_text(encoding="utf-8", errors="replace").splitlines()[0].strip()
            except Exception:
                first_line = None
            entry["agents_first_line"] = first_line

        # Policy suggestion (operator-readable; not a hard law).
        if name == "trench_re_vault":
            entry["suggested_role"] = "DIRTY airlock (reverse-engineering / capture / extraction)."
            entry["suggested_policy"] = "Read-only from CLEAN/shipping code. Promote only sanitized contracts into this repo."
        elif name == "trenchwork_clean":
            entry["suggested_role"] = "Separate CLEAN workspace snapshot (factory/tooling)."
            entry["suggested_policy"] = "Treat as non-canonical here; copy-in only specific files with explicit intent + provenance."
        else:
            entry["suggested_role"] = "external workspace"
            entry["suggested_policy"] = "non-canonical"

        out.append(entry)

    return out


def _write_truth_json(path: pathlib.Path) -> dict[str, Any]:
    shape = _shape_counts()
    p2k = _p2k_counts()
    external = _external_workspaces()
    truth_files = _scan_truth_files()

    truth: dict[str, Any] = {
        "schema": "trench.repo_truth.v1",
        "generated_at": _iso_now(),
        "repo_root_abs": str(REPO_ROOT),
        "canonical": {
            "agents": "AGENTS.md",
            "claude": "CLAUDE.md",
            "spec": "SPEC.md",
            "shipping": "SHIPPING.md",
            "repo_truth_md": "REPO_TRUTH.md",
            "codebase_map": "CODEBASE_MAP.md",
            "dsp_rules": "trench-core/CLAUDE.md",
            "pyruntime_rules": "pyruntime/CLAUDE.md",
            "shapes_readme": "vault/_shapes/README.md",
            "phonemes_readme": "vault/_phonemes/README.md",
        },
        "entrypoints": {
            "rust_workspace": "Cargo.toml",
            "python_api": "pyruntime/api.py",
            "python_tools_dir": "tools/",
            "juce_repo_dir": "trench-juce/",
            "juce_plugin_repo_dir": "trench-juce/plugin/",
        },
        "verification": {
            "rust": [
                "cargo test --workspace",
                "cargo clippy --workspace --all-targets -- -D warnings",
            ],
            "python": [
                "python -m pytest tests/test_api.py -v",
            ],
            "shape_bank": [
                "python tools/author_sonic_bank.py --pairs --out-root vault/_shapes",
            ],
        },
        "data_products": {
            "sonic_tables": {
                "tables": "docs/sonic_tables/tables.json",
                "tables_bark": "docs/sonic_tables/tables_bark.json",
                "heritage_phonemes": "docs/sonic_tables/heritage_phonemes.json",
                "unified_bark": "docs/sonic_tables/sonic_tables_unified_bark.json",
            },
            "shape_bank": {
                "root": "vault/_shapes/",
                "manifest": "vault/_shapes/manifest.json",
                "manifest_unified": "vault/_shapes/manifest_unified.json",
                "catalog": "vault/_shapes/shape_bank_catalog.json",
                "catalog_unified": "vault/_shapes/shape_bank_catalog_unified.json",
                "palette": "vault/_shapes/palette.jsonl",
                "palette_unified": "vault/_shapes/palette_unified.jsonl",
                "counts": shape.get("counts", {}),
                "counts_source": shape.get("source"),
            },
            "phonemes": {
                "p2k_inventory_raw": "vault/_phonemes/p2k_phoneme_inventory.json",
                "p2k_inventory_v2": "vault/_phonemes/p2k_phoneme_inventory_v2.json",
                "p2k_aliases_v2": "vault/_phonemes/p2k_phoneme_aliases_v2.json",
                "token_inventory_unified_v2": "vault/_phonemes/token_inventory_unified_v2.json",
                "unified_sonic_inventory_v2": "vault/_phonemes/unified_sonic_inventory_v2.json",
                "counts": p2k,
            },
        },
        "external_workspaces": external,
        "truth_files": truth_files,
        "non_canonical_dirs": [
            "_GRAVEYARD/",
            "dev/tmp/",
            "target/",
            "target-codex-vizia/",
        ],
        "existence_checks": {
            # Cheap sanity: if these flip to false, the index is stale.
            "SPEC.md": _exists("SPEC.md"),
            "SHIPPING.md": _exists("SHIPPING.md"),
            "Cargo.toml": _exists("Cargo.toml"),
            "pyruntime/api.py": _exists("pyruntime/api.py"),
            "trench-juce/plugin": _exists("trench-juce/plugin"),
            "vault/_shapes/manifest.json": _exists("vault/_shapes/manifest.json"),
            "vault/_phonemes/p2k_phoneme_inventory_v2.json": _exists("vault/_phonemes/p2k_phoneme_inventory_v2.json"),
        },
        "notes": [
            "SPEC.md is the immutable law. REPO_TRUTH.json is the live index; update via tools/update_repo_truth.py after moving/renaming paths.",
            "If docs disagree with build/runtime wiring, report the conflict explicitly and follow the build-wired truth.",
        ],
    }

    path.write_text(json.dumps(truth, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return truth


def _write_truth_md(path: pathlib.Path, truth: dict[str, Any]) -> None:
    shape = (truth.get("data_products", {}).get("shape_bank", {}) or {})
    shape_counts = shape.get("counts") or {}
    p2k = (truth.get("data_products", {}).get("phonemes", {}).get("counts") or {})

    def _line(s: str) -> str:
        return s.rstrip() + "\n"

    lines: list[str] = []
    lines.append(_line("# Repo Truth (Canonical Index)"))
    lines.append(_line(""))
    lines.append(_line("This is the single jump-table for **what is canonical** in the repo."))
    lines.append(_line("If paths/entrypoints change, regenerate with: `python tools/update_repo_truth.py`."))
    lines.append(_line(""))

    lines.append(_line("## Start Here"))
    lines.append(_line("- Product/DSP law: `SPEC.md`"))
    lines.append(_line("- Shipping checklist + release gate: `SHIPPING.md`"))
    lines.append(_line("- Agent posture + bans: `AGENTS.md`, `CLAUDE.md`"))
    lines.append(_line("- DSP subtree rules: `trench-core/CLAUDE.md`"))
    lines.append(_line("- Forge/compiler boundary rules: `pyruntime/CLAUDE.md`"))
    lines.append(_line(""))

    lines.append(_line("## Entrypoints"))
    ep = truth.get("entrypoints", {}) or {}
    for k in ("rust_workspace", "python_api", "juce_plugin_repo_dir"):
        if k in ep:
            lines.append(_line(f"- {k}: `{ep[k]}`"))
    lines.append(_line(""))

    lines.append(_line("## Verification (Actual Commands)"))
    ver = truth.get("verification", {}) or {}
    if isinstance(ver.get("rust"), list):
        lines.append(_line("- Rust:"))
        for cmd in ver["rust"]:
            lines.append(_line(f"  - `{cmd}`"))
    if isinstance(ver.get("python"), list):
        lines.append(_line("- Python:"))
        for cmd in ver["python"]:
            lines.append(_line(f"  - `{cmd}`"))
    if isinstance(ver.get("shape_bank"), list):
        lines.append(_line("- Shape bank regen:"))
        for cmd in ver["shape_bank"]:
            lines.append(_line(f"  - `{cmd}`"))
    lines.append(_line(""))

    lines.append(_line("## Shape Bank (Current)"))
    lines.append(_line(f"- counts_source: `{shape.get('counts_source')}`"))
    if isinstance(shape_counts, dict) and shape_counts:
        lines.append(_line(f"- total_bodies: {shape_counts.get('total_bodies')}"))
        lines.append(_line(f"- static_bodies: {shape_counts.get('static_bodies')}"))
        lines.append(_line(f"- pair_bodies: {shape_counts.get('pair_bodies')}"))
        by_cat = shape_counts.get("by_category")
        if isinstance(by_cat, dict):
            lines.append(_line(f"- categories: {len(by_cat)}"))
    lines.append(_line("- bank docs: `vault/_shapes/README.md`"))
    lines.append(_line(""))

    lines.append(_line("## P2K Phonemes (Current)"))
    lines.append(_line(f"- clusters: {p2k.get('clusters')}"))
    lines.append(_line(f"- unmapped: {p2k.get('unmapped')}"))
    lines.append(_line(f"- heritage_fills: {p2k.get('heritage_fills')}"))
    lines.append(_line("- inventory docs: `vault/_phonemes/README.md`"))
    lines.append(_line(""))

    tf = truth.get("truth_files", {}) or {}
    if isinstance(tf, dict) and tf:
        lines.append(_line("## Truth Files (Single Voice)"))
        lines.append(_line(f"- ok: {tf.get('ok')}"))
        v = tf.get("violations")
        if isinstance(v, list):
            lines.append(_line(f"- violations: {len(v)}"))
        found = tf.get("found")
        rules = tf.get("rules") or {}
        if isinstance(found, dict) and isinstance(rules, dict):
            for name in ("SPEC.md", "SHIPPING.md", "AGENTS.md", "CLAUDE.md"):
                if name not in found:
                    continue
                canonical = (rules.get(name) or {}).get("canonical")
                lines.append(_line(f"- {name}: canonical `{canonical}`"))
                # List stubs/allowed if present.
                for e in found.get(name) or []:
                    if not isinstance(e, dict):
                        continue
                    kind = e.get("kind")
                    if kind in ("canonical",):
                        continue
                    p = e.get("path")
                    ok = e.get("ok")
                    lines.append(_line(f"  - {kind} ok={ok}: `{p}`"))
        lines.append(_line(""))

    lines.append(_line("## Non-canonical (Don’t Treat As Truth)"))
    for d in truth.get("non_canonical_dirs", []) or []:
        lines.append(_line(f"- `{d}`"))
    lines.append(_line(""))

    ext = truth.get("external_workspaces", []) or []
    if isinstance(ext, list) and ext:
        lines.append(_line("## External Workspaces (Non-canonical)"))
        lines.append(
            _line(
                "These directories are outside this repo. They may contain related work, but they are not the shipping truth for `C:\\Users\\hooki\\Trench`."
            )
        )
        for e in ext:
            if not isinstance(e, dict):
                continue
            name = e.get("name")
            path_abs = e.get("path_abs")
            exists = e.get("exists")
            role = e.get("suggested_role")
            policy = e.get("suggested_policy")
            if name and path_abs:
                lines.append(_line(f"- `{name}`: `{path_abs}` (exists={exists})"))
                if role:
                    lines.append(_line(f"  - role: {role}"))
                if policy:
                    lines.append(_line(f"  - policy: {policy}"))
        lines.append(_line(""))

    lines.append(_line("## Existence Checks (Staleness Alarm)"))
    checks = truth.get("existence_checks", {}) or {}
    if isinstance(checks, dict):
        for k in sorted(checks.keys()):
            lines.append(_line(f"- `{k}`: {checks[k]}"))

    path.write_text("".join(lines), encoding="utf-8")


REPO_ROOT = pathlib.Path(__file__).resolve().parent.parent


def main() -> int:
    truth_json_path = REPO_ROOT / "REPO_TRUTH.json"
    truth_md_path = REPO_ROOT / "REPO_TRUTH.md"

    truth = _write_truth_json(truth_json_path)
    _write_truth_md(truth_md_path, truth)

    print(f"WROTE {truth_json_path}")
    print(f"WROTE {truth_md_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
