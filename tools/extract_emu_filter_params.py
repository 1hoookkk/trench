"""Extract E-mu MorphDesigner filter parameters from NotebookLM markdown exports.

The NotebookLM bundle (`NotebookLM_Trench_Bundle.md`, gitignored at
.gitignore:23) and its split files `10-templates-filter-01.md` …
`13-templates-filter-04.md` ship E-mu's authoring XML in plaintext bullet
form, e.g.::

    - hedz/filter/type-absolute: 144
    - hedz/filter/designer-section[1]/type: 3
    - hedz/filter/designer-section[1]/low-freq: 0
    - hedz/filter/designer-section[1]/low-gain: 127
    - hedz/filter/designer-section[1]/high-freq: 116
    - hedz/filter/designer-section[1]/high-gain: 87

That is the exact 30-integer authoring grid (6 stages × 5 fields) the
COMPILER side of `pyruntime/designer_compile.py` already consumes via
`parse_xml`. This tool is the markdown sibling of `parse_xml`: it walks
a path of NotebookLM `.md` exports, groups bullets by template name +
section index, and writes a JSON inventory that
`pyruntime.designer_compile.load_inventory` reads back into
`DesignerTemplate` records.

CLI:

    python tools/extract_emu_filter_params.py \\
        --input <path-to-md-or-dir> \\
        --out vault/_phonemes/heritage_designer_sections.json

Stdlib only — no pyruntime imports, no third-party deps. Stays runnable
from a clean venv on the user's Windows machine where the NotebookLM
bundle actually lives.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

INVENTORY_FORMAT = "heritage-designer-sections-v1"

# One authoritative bullet pattern. Captures:
#   name  — template identifier (anything up to the first '/')
#   rest  — everything between '/filter/' and ':'
#   value — integer or float right-hand side
_BULLET_RE = re.compile(
    r"^\s*[-*]?\s*(?P<name>[A-Za-z0-9_][^/\s]*)/filter/(?P<rest>\S+?):\s*(?P<value>-?\d+(?:\.\d+)?)\s*$"
)

_SECTION_RE = re.compile(
    r"^designer-section\[(?P<idx>[1-6])\]/(?P<field>type|low-freq|low-gain|high-freq|high-gain)$"
)

_FIELD_KEY = {
    "type": "type",
    "low-freq": "low_freq",
    "low-gain": "low_gain",
    "high-freq": "high_freq",
    "high-gain": "high_gain",
}


def _empty_section(idx: int) -> dict:
    return {
        "index": idx,
        "type": 0,
        "low_freq": 0,
        "low_gain": 0,
        "high_freq": 0,
        "high_gain": 0,
    }


@dataclass
class RawTemplate:
    """In-memory accumulator for one NotebookLM filter template.

    Mirrors `pyruntime.designer_compile.DesignerTemplate` semantics
    (1-indexed sections, 6 sections per template, missing sections
    zero-fill) without importing pyruntime — keeps this tool runnable
    from a clean venv.
    """

    name: str
    type_absolute: int | None = None
    frequency: float = 0.0
    gain: float = 0.0
    sections: list[dict] = field(
        default_factory=lambda: [_empty_section(i) for i in range(1, 7)]
    )

    def set_section_field(self, idx: int, field_name: str, value: int) -> None:
        # idx is 1..6, list is 0..5
        self.sections[idx - 1][_FIELD_KEY[field_name]] = int(value)


def parse_notebooklm_export(text: str) -> "OrderedDict[str, RawTemplate]":
    """Parse a single NotebookLM markdown blob into RawTemplate records.

    Templates are returned in first-seen order so the inventory is stable
    across reruns of the same input. Bullets that don't match the
    `<name>/filter/...` shape are silently skipped — NotebookLM exports
    interleave prose, headings, and other bullets we don't care about.
    """
    templates: "OrderedDict[str, RawTemplate]" = OrderedDict()

    for line in text.splitlines():
        m = _BULLET_RE.match(line)
        if not m:
            continue

        name = m.group("name")
        rest = m.group("rest")
        raw_value = m.group("value")

        tpl = templates.get(name)
        if tpl is None:
            tpl = RawTemplate(name=name)
            templates[name] = tpl

        sec_match = _SECTION_RE.match(rest)
        if sec_match:
            idx = int(sec_match.group("idx"))
            field_name = sec_match.group("field")
            tpl.set_section_field(idx, field_name, int(float(raw_value)))
            continue

        if rest == "type-absolute":
            tpl.type_absolute = int(float(raw_value))
        elif rest == "frequency":
            tpl.frequency = float(raw_value)
        elif rest == "gain":
            tpl.gain = float(raw_value)
        # else: unknown filter sub-key, ignore

    return templates


def extract_directory(root: Path) -> "OrderedDict[str, RawTemplate]":
    """Walk a path (file or directory) of NotebookLM .md exports.

    Files are visited in sorted order so multi-file extractions are
    deterministic. The first appearance of a template name wins; later
    appearances merge their fields onto the existing record so split
    bundles still recombine cleanly.
    """
    templates: "OrderedDict[str, RawTemplate]" = OrderedDict()

    if root.is_file():
        files: Iterable[Path] = [root]
    else:
        files = sorted(p for p in root.rglob("*.md") if p.is_file())

    for path in files:
        text = path.read_text(encoding="utf-8", errors="replace")
        for name, tpl in parse_notebooklm_export(text).items():
            existing = templates.get(name)
            if existing is None:
                templates[name] = tpl
                continue
            # Merge: later wins on populated fields only
            if tpl.type_absolute is not None:
                existing.type_absolute = tpl.type_absolute
            if tpl.frequency:
                existing.frequency = tpl.frequency
            if tpl.gain:
                existing.gain = tpl.gain
            for i in range(6):
                src = tpl.sections[i]
                # Only overwrite the section if the incoming one is non-empty
                if any(src[k] for k in ("type", "low_freq", "low_gain", "high_freq", "high_gain")):
                    existing.sections[i] = src

    return templates


def shape_id_for(template: RawTemplate) -> str:
    """Stable shape_id matching the convention in
    docs/looperator/filter_cube_display_model.md:62 (`abs<n>_<hex12>`).

    Hashes the 30-integer section tuple only — frequency/gain defaults
    are display-time parameters per the same doc and must not perturb
    the structural identity.
    """
    payload_tuple = tuple(
        (s["type"], s["low_freq"], s["low_gain"], s["high_freq"], s["high_gain"])
        for s in template.sections
    )
    digest = hashlib.sha256(repr(payload_tuple).encode("utf-8")).hexdigest()[:12]
    abs_key = template.type_absolute if template.type_absolute is not None else 0
    return f"abs{abs_key}_{digest}"


def to_inventory(
    templates: "OrderedDict[str, RawTemplate]",
    extracted_from: str = "",
) -> dict:
    """Render templates into the heritage-designer-sections-v1 schema."""
    out_templates = []
    for tpl in templates.values():
        out_templates.append(
            {
                "name": tpl.name,
                "type_absolute": tpl.type_absolute,
                "frequency": tpl.frequency,
                "gain": tpl.gain,
                "shape_id": shape_id_for(tpl),
                "sections": [dict(s) for s in tpl.sections],
            }
        )
    return {
        "format": INVENTORY_FORMAT,
        "extracted_from": extracted_from,
        "template_count": len(out_templates),
        "templates": out_templates,
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Extract E-mu MorphDesigner filter parameters from NotebookLM "
            "markdown exports into a JSON inventory consumable by "
            "pyruntime.designer_compile.load_inventory."
        ),
    )
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="Path to a NotebookLM .md file or a directory containing .md exports.",
    )
    parser.add_argument(
        "--out",
        required=True,
        type=Path,
        help="Output JSON inventory path.",
    )
    parser.add_argument(
        "--per-template",
        type=Path,
        default=None,
        help=(
            "Optional directory; when set, each template is also written as "
            "a standalone JSON file <name>.json alongside the combined inventory."
        ),
    )
    args = parser.parse_args(argv)

    if not args.input.exists():
        parser.error(f"--input does not exist: {args.input}")

    templates = extract_directory(args.input)
    if not templates:
        parser.error(f"no NotebookLM filter bullets found under {args.input}")

    inventory = to_inventory(templates, extracted_from=str(args.input))

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(inventory, indent=2, sort_keys=False) + "\n")

    if args.per_template is not None:
        args.per_template.mkdir(parents=True, exist_ok=True)
        for entry in inventory["templates"]:
            (args.per_template / f"{entry['name']}.json").write_text(
                json.dumps(entry, indent=2, sort_keys=False) + "\n"
            )

    print(f"extracted {inventory['template_count']} templates → {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
