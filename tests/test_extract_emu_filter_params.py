"""Tests for tools/extract_emu_filter_params.py.

Stdlib unittest only — no pytest dep — so this runs in `./check`
without expanding the toolchain. Run directly with::

    python -m unittest tests.test_extract_emu_filter_params -v

The fixtures under `tests/fixtures/notebooklm/` are the real NotebookLM
markdown exports of the E-mu MorphDesigner XML (10-13-templates-filter-
*.md, 83 templates total). The Talking Hedz stage-1 anchor and the
three documented `shape_id` alias groups from
`docs/looperator/filter_cube_display_model.md:60-64` are the ground
truth this test pins against. If any anchor or alias-group hash
changes, the parser or the canonical hash scheme has drifted and
needs reconciliation, not "tolerance".
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
from collections import defaultdict
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tools.extract_emu_filter_params import (  # noqa: E402
    INVENTORY_FORMAT,
    RawTemplate,
    extract_directory,
    parse_notebooklm_export,
    shape_id_for,
    to_inventory,
)

# pyruntime pulls in numpy via pyruntime.encode. The full ./check
# environment has numpy; bare sandboxes may not. Skip the
# designer_compile round-trip if pyruntime can't be imported, so the
# parser-only tests stay green everywhere.
try:
    from pyruntime import designer_compile as _designer_compile  # noqa: F401
    _PYRUNTIME_AVAILABLE = True
except Exception:
    _PYRUNTIME_AVAILABLE = False

NOTEBOOKLM_DIR = REPO_ROOT / "tests" / "fixtures" / "notebooklm"
HEDZ_FILE = NOTEBOOKLM_DIR / "13-templates-filter-04.md"

# Documented canonical alias groups — the only validation anchor for
# the shape_id hash scheme. Source: docs/looperator/filter_cube_display_model.md:60-64.
EXPECTED_ALIAS_GROUPS = {
    "abs108_b352f374685b": {"AC Wow", "AC Morph"},
    "abs144_b5c825095c87": {"50", "100", "hedz0", "d", "hedz"},
    "abs144_7949e309f941": {"morph0q100", "morph100q100", "ssssss"},
}


class ParseNotebookLMExportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.text = HEDZ_FILE.read_text(encoding="utf-8")
        self.templates = parse_notebooklm_export(self.text)

    def test_hedz_stage1_locks_user_quoted_truth(self) -> None:
        """Stage 1 of `hedz` must equal the user-quoted Talking Hedz
        `13-templates-filter-04.md` values verbatim. If this fails the
        parser regressed against the only ground-truth anchor we have
        from the user's original task description."""
        tpl = self.templates["hedz"]
        s1 = tpl.sections[0]
        self.assertEqual(s1["index"], 1)
        self.assertEqual(s1["type"], 3)
        self.assertEqual(s1["low_freq"], 0)
        self.assertEqual(s1["low_gain"], 127)
        self.assertEqual(s1["high_freq"], 116)
        self.assertEqual(s1["high_gain"], 87)

    def test_type_absolute_and_filter_defaults(self) -> None:
        hedz = self.templates["hedz"]
        self.assertEqual(hedz.type_absolute, 144)
        self.assertAlmostEqual(hedz.frequency, 1.0)
        self.assertAlmostEqual(hedz.gain, 1.0)

    def test_all_six_sections_present(self) -> None:
        hedz = self.templates["hedz"]
        self.assertEqual(len(hedz.sections), 6)
        for i, sec in enumerate(hedz.sections, start=1):
            self.assertEqual(sec["index"], i)

    def test_multi_word_template_names_parse(self) -> None:
        """Real E-mu names contain spaces and digits — the parser
        must not split on them. File 13 contains `morph0q100`,
        `morph100q100`, `ssssss`. File 10 contains `All Off Designer`
        and `Bass Boost LP 1`. Verify both shapes parse from their
        respective files."""
        # File 13 templates
        for name in ("hedz", "hedz0", "morph0q100", "morph100q100", "ssssss"):
            self.assertIn(name, self.templates, f"missing {name!r} from file 13")

        text10 = (NOTEBOOKLM_DIR / "10-templates-filter-01.md").read_text(encoding="utf-8")
        templates10 = parse_notebooklm_export(text10)
        for name in ("All Off Designer", "Bass Boost LP 1"):
            self.assertIn(name, templates10, f"missing {name!r} from file 10")

    def test_unrelated_lines_are_ignored(self) -> None:
        """Prose, headings, `@type: long` schema bullets, and bullets
        that don't match the `<name>/filter/...:<number>` shape must
        be silently dropped."""
        noisy = (
            "# Some heading\n"
            "Some prose paragraph.\n"
            "- not_a_filter_bullet: 42\n"
            "- foo/notfilter/something: 5\n"
            "- bar/filter/designer-section[1]/type @type: long\n"
            "- bar/filter/designer-section[1]/type: 7\n"
        )
        templates = parse_notebooklm_export(noisy)
        self.assertEqual(list(templates.keys()), ["bar"])
        self.assertEqual(templates["bar"].sections[0]["type"], 7)


class ShapeIdTests(unittest.TestCase):
    def _make_template(self) -> RawTemplate:
        tpl = RawTemplate(name="x", type_absolute=144)
        tpl.set_section_field(1, "type", 3)
        tpl.set_section_field(1, "low-freq", 0)
        tpl.set_section_field(1, "low-gain", 127)
        tpl.set_section_field(1, "high-freq", 116)
        tpl.set_section_field(1, "high-gain", 87)
        return tpl

    def test_shape_id_is_deterministic(self) -> None:
        a = self._make_template()
        b = self._make_template()
        self.assertEqual(shape_id_for(a), shape_id_for(b))

    def test_shape_id_namespaces_by_type_absolute(self) -> None:
        a = self._make_template()
        b = self._make_template()
        b.type_absolute = 108
        self.assertNotEqual(shape_id_for(a), shape_id_for(b))
        self.assertTrue(shape_id_for(a).startswith("abs144_"))
        self.assertTrue(shape_id_for(b).startswith("abs108_"))

    def test_shape_id_ignores_frequency_and_gain_defaults(self) -> None:
        """Per docs/looperator/filter_cube_display_model.md:69
        frequency/gain are display-time defaults, not structural
        identity. Two templates that differ only in those defaults
        must share a shape_id."""
        a = self._make_template()
        b = self._make_template()
        b.frequency = 0.99
        b.gain = -0.5
        self.assertEqual(shape_id_for(a), shape_id_for(b))


class InventoryTests(unittest.TestCase):
    def setUp(self) -> None:
        self.templates = extract_directory(NOTEBOOKLM_DIR)
        self.inventory = to_inventory(self.templates, extracted_from="fixtures")

    def test_inventory_format_and_count(self) -> None:
        self.assertEqual(self.inventory["format"], INVENTORY_FORMAT)
        # 83 templates in 10-13-templates-filter-*.md (file headers say
        # "Templates in this file: 25/25/25/8 of 83").
        self.assertEqual(self.inventory["template_count"], 83)
        self.assertEqual(self.inventory["extracted_from"], "fixtures")

    def test_inventory_round_trips_through_json(self) -> None:
        blob = json.dumps(self.inventory)
        decoded = json.loads(blob)
        self.assertEqual(decoded["format"], INVENTORY_FORMAT)
        hedz = next(t for t in decoded["templates"] if t["name"] == "hedz")
        self.assertEqual(hedz["sections"][0]["type"], 3)
        self.assertEqual(hedz["sections"][0]["low_gain"], 127)
        self.assertEqual(hedz["sections"][0]["high_freq"], 116)
        self.assertEqual(hedz["sections"][0]["high_gain"], 87)

    def test_canonical_alias_groups_reproduce(self) -> None:
        """The three documented `shape_id` alias groups from
        docs/looperator/filter_cube_display_model.md:60-64 must
        reproduce exactly. If any of these break, the canonical hash
        scheme has drifted — fix `shape_id_for`, do not relax the test."""
        groups: dict[str, set[str]] = defaultdict(set)
        for entry in self.inventory["templates"]:
            groups[entry["shape_id"]].add(entry["name"])

        for shape_id, expected_names in EXPECTED_ALIAS_GROUPS.items():
            actual = groups.get(shape_id, set())
            self.assertEqual(
                actual,
                expected_names,
                f"alias group {shape_id} drifted: expected {sorted(expected_names)} "
                f"got {sorted(actual)}",
            )

    def test_type_absolute_distribution_matches_source(self) -> None:
        """Files 10-12 are 25 templates each, file 13 is 8 → 83 total.
        Distribution against the real markdown should be 108=72,
        144=9, 140=2. The 4-template gap vs the canonical
        analysis script's count of 79 (108=70, 144=9) is the four
        fully-zeroed bypass templates the canonical filtered out:
        `All Off Designer`, `des` (108) and `tb`, `w` (140). The
        parser intentionally keeps them so the inventory is faithful
        to the source markdown."""
        ta_counts: dict[int | None, int] = defaultdict(int)
        for entry in self.inventory["templates"]:
            ta_counts[entry["type_absolute"]] += 1
        self.assertEqual(ta_counts[108], 72)
        self.assertEqual(ta_counts[144], 9)
        self.assertEqual(ta_counts[140], 2)

    @unittest.skipUnless(_PYRUNTIME_AVAILABLE, "pyruntime (numpy) not importable")
    def test_designer_compile_load_inventory_round_trip(self) -> None:
        """End-to-end: inventory JSON → designer_compile.load_inventory
        → DesignerTemplate → compile_designer must succeed for every
        template in the fixture without raising."""
        from pyruntime import designer_compile

        with tempfile.TemporaryDirectory() as tmp:
            inv_path = Path(tmp) / "inv.json"
            inv_path.write_text(json.dumps(self.inventory))

            templates = designer_compile.load_inventory(inv_path)
            self.assertEqual(len(templates), 83)

            for tpl in templates:
                self.assertEqual(len(tpl.sections), 6)
                # Every section must be a real DesignerSection with the
                # five integer fields populated (zero is a valid value).
                for sec in tpl.sections:
                    self.assertIsInstance(sec.type, int)
                    self.assertIsInstance(sec.low_freq, int)
                    self.assertIsInstance(sec.low_gain, int)
                    self.assertIsInstance(sec.high_freq, int)
                    self.assertIsInstance(sec.high_gain, int)

                corners = designer_compile.compile_designer(tpl)
                # Four corners (M0/M100 × Q0/Q100) by CornerArray contract.
                self.assertIsNotNone(corners.a)
                self.assertIsNotNone(corners.b)
                self.assertIsNotNone(corners.c)
                self.assertIsNotNone(corners.d)


class ExtractDirectoryTests(unittest.TestCase):
    def test_single_file_extraction(self) -> None:
        templates = extract_directory(HEDZ_FILE)
        self.assertIn("hedz", templates)
        # File 13 header: "Templates in this file: 8 of 83".
        self.assertEqual(len(templates), 8)

    def test_directory_walk_extracts_all_83_templates(self) -> None:
        templates = extract_directory(NOTEBOOKLM_DIR)
        self.assertEqual(len(templates), 83)

    def test_directory_walk_merges_split_bundles(self) -> None:
        """Templates split across files (the NotebookLM 10-13 case)
        must merge cleanly: the second file's populated fields land on
        the first file's record."""
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp)
            (d / "01.md").write_text(
                "- split_target/filter/type-absolute: 144\n"
                "- split_target/filter/designer-section[1]/type: 3\n"
                "- split_target/filter/designer-section[1]/low-freq: 0\n"
                "- split_target/filter/designer-section[1]/low-gain: 127\n"
                "- split_target/filter/designer-section[1]/high-freq: 116\n"
                "- split_target/filter/designer-section[1]/high-gain: 87\n"
            )
            (d / "02.md").write_text(
                "- split_target/filter/designer-section[6]/type: 1\n"
                "- split_target/filter/designer-section[6]/low-freq: 60\n"
                "- split_target/filter/designer-section[6]/low-gain: 60\n"
                "- split_target/filter/designer-section[6]/high-freq: 66\n"
                "- split_target/filter/designer-section[6]/high-gain: 66\n"
            )
            templates = extract_directory(d)

        self.assertEqual(list(templates.keys()), ["split_target"])
        tpl = templates["split_target"]
        self.assertEqual(tpl.type_absolute, 144)
        self.assertEqual(tpl.sections[0]["type"], 3)
        self.assertEqual(tpl.sections[0]["high_gain"], 87)
        self.assertEqual(tpl.sections[5]["type"], 1)
        self.assertEqual(tpl.sections[5]["high_freq"], 66)


if __name__ == "__main__":
    unittest.main()
