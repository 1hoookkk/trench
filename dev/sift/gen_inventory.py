"""Build script_inventory.md for Trench, trench_re_vault, trenchwork_clean."""
import os, re, subprocess

ROOTS = [
    "C:/Users/hooki/Trench",
    "C:/Users/hooki/trench_re_vault",
    "C:/Users/hooki/trenchwork_clean",
]

SKIP = {
    "/build/", "/build-codex/", "/Builds/", "/JUCE/", "/_deps/",
    "/_graveyard/", "/.claude/worktrees/", "/catch2-src/",
}

def skip(path):
    p = path.replace("\\", "/")
    return any(s in p for s in SKIP)

def find_scripts(root):
    py, rs = [], []
    for dirpath, dirs, files in os.walk(root):
        dp = dirpath.replace("\\", "/")
        if any(s in dp for s in SKIP):
            dirs[:] = []
            continue
        for f in files:
            full = os.path.join(dirpath, f).replace("\\", "/")
            if f.endswith(".py"):
                py.append(full)
            elif f == "main.rs" or (f.endswith(".rs") and "/src/bin/" in full):
                rs.append(full)
    # Filter py: must have main/__main__/argparse
    exe_py = []
    for p in py:
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as fh:
                txt = fh.read(3000)
            if re.search(r'def main\(|__main__|argparse', txt):
                exe_py.append(p)
        except Exception:
            pass
    return exe_py, rs

def get_purpose(path, rust=False):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            txt = f.read(2000)
        if not rust:
            for pat in [r'"""(.+?)"""', r"'''(.+?)'''"]:
                m = re.search(pat, txt, re.DOTALL)
                if m:
                    s = m.group(1).strip().split("\n")[0].strip()
                    if len(s) > 5:
                        return s[:100]
        for line in txt.split("\n")[:20]:
            line = line.strip()
            prefix = "//" if rust else "#"
            if line.startswith(prefix) and not line.startswith(prefix * 3) and not line.startswith("#!"):
                c = line.lstrip(prefix).strip()
                if len(c) > 5:
                    return c[:100]
        return "(no docstring)"
    except Exception:
        return "(unreadable)"

def git_date(path):
    d, b = os.path.dirname(path), os.path.basename(path)
    try:
        r = subprocess.run(["git", "-C", d, "log", "-1", "--format=%as", "--", b],
                           capture_output=True, text=True, timeout=5)
        s = r.stdout.strip()
        if s:
            return s
    except Exception:
        pass
    try:
        import datetime
        return datetime.date.fromtimestamp(os.path.getmtime(path)).isoformat()
    except Exception:
        return "?"

def writes_files(path):
    try:
        with open(path, "r", encoding="utf-8", errors="ignore") as f:
            txt = f.read()
        if re.search(r'open\([^)]+["\']w["\']|\.write_text\(|json\.dump\(|Path\([^)]+\)\.write', txt):
            outs = re.findall(r'["\']([^"\']{5,60}\.(?:json|yaml|bin|txt|csv|png|md))["\']', txt)
            outs = [x for x in outs if not x.startswith("http")][:2]
            return "yes — " + ", ".join(outs) if outs else "yes"
        return "no"
    except Exception:
        return "?"

def is_referenced(basename, all_py, all_rs):
    stem = os.path.splitext(basename)[0]
    for p in all_py + all_rs:
        if os.path.basename(p) == basename:
            continue
        try:
            with open(p, "r", encoding="utf-8", errors="ignore") as f:
                if stem in f.read():
                    rel = os.path.relpath(p, "C:/Users/hooki")
                    return f"yes — {rel}"
        except Exception:
            pass
    return "no"

HEADER = "| Path | Purpose | Last Modified | Referenced By | Writes Files |\n|------|---------|---------------|---------------|--------------|"

# Collect all files first so cross-ref check works
all_py_global, all_rs_global = [], []
per_root = {}
for root in ROOTS:
    py, rs = find_scripts(root)
    per_root[root] = (py, rs)
    all_py_global += py
    all_rs_global += rs

# Build output
lines = ["# Script Inventory\n", "_Generated 2026-04-17_\n"]

all_py_for_ref = all_py_global  # use py files only for ref check (fast enough)

for root in ROOTS:
    py, rs = per_root[root]
    label = root.replace("C:/Users/hooki/", "")
    lines.append(f"\n## {root}\n")

    if py:
        lines.append(f"### Python ({len(py)} scripts)\n")
        lines.append(HEADER)
        for p in sorted(py):
            purpose = get_purpose(p)
            date = git_date(p)
            ref = is_referenced(os.path.basename(p), all_py_for_ref, all_rs_global)
            writes = writes_files(p)
            lines.append(f"| `{p}` | {purpose} | {date} | {ref} | {writes} |")

    if rs:
        lines.append(f"\n### Rust binaries ({len(rs)})\n")
        lines.append(HEADER)
        for p in sorted(rs):
            purpose = get_purpose(p, rust=True)
            date = git_date(p)
            ref = is_referenced(os.path.basename(p), all_py_for_ref, all_rs_global)
            lines.append(f"| `{p}` | {purpose} | {date} | {ref} | — |")

# Duplicates
lines.append("\n## Duplicate stems\n")
stems = {}
for p in all_py_global:
    s = os.path.splitext(os.path.basename(p))[0]
    stems.setdefault(s, []).append(p)
dupes = {k: v for k, v in stems.items() if len(v) > 1}
if dupes:
    for stem, paths in sorted(dupes.items()):
        lines.append(f"- **{stem}** ({len(paths)}x): " + " | ".join(paths))
else:
    lines.append("None found.")

out_path = "C:/Users/hooki/Trench/dev/sift/script_inventory.md"
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w", encoding="utf-8") as f:
    f.write("\n".join(lines))

total_py = sum(len(v[0]) for v in per_root.values())
total_rs = sum(len(v[1]) for v in per_root.values())
print(f"Done: {total_py} Python + {total_rs} Rust = {total_py+total_rs} total")
print(f"Duplicate stems: {sorted(dupes.keys())}")
