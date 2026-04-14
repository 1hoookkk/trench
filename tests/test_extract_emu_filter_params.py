"""Tests for tools/extract_emu_filter_params.py.

Stdlib unittest only — no pytest dep — so this runs in `./check`
without expanding the toolchain. Run directly with::

    python -m unittest tests.test_extract_emu_filter_params -v
"""
from __future__ import annotations

import json
import sys
import tempfile
import unittest
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

FIXTURE = REPO_ROOT / "tests" / "fixtures" / "notebooklm_hedz_excerpt.md"


class ParseNotebookLMExportTests(unittest.TestCase):
    def setUp(self) -> None:
        self.text = FIXTURE.read_text(encoding="utf-8")
        self.templates = parse_notebooklm_export(self.text)

    def test_hedz_stage1_locks_user_quoted_truth(self) -> None:
        """Stage 1 of synthetic_hedz must equal the user-quoted Talking
        Hedz `13-templates-filter-04.md` values verbatim. If this fails
        the parser broke the only real anchor in the fixture."""
        tpl = self.templates["synthetic_hedz"]
        s1 = tpl.sections[0]
        self.assertEqual(s1["index"], 1)
        self.assertEqual(s1["type"], 3)
        self.assertEqual(s1["low_freq"], 0)
        self.assertEqual(s1["low_gain"], 127)
        self.assertEqual(s1["high_freq"], 116)
        self.assertEqual(s1["high_gain"], 87)

    def test_type_absolute_and_filter_defaults(self) -> None:
        hedz = self.templates["synthetic_hedz"]
        self.assertEqual(hedz.type_absolute, 144)
        self.assertEqual(hedz.frequency, 0.0)
        self.assertEqual(hedz.gain, 0.0)

        acwow = self.templates["synthetic_acwow"]
        self.assertEqual(acwow.type_absolute, 108)
        self.assertAlmostEqual(acwow.frequency, 0.5)
        self.assertAlmostEqual(acwow.gain, 0.25)

    def test_all_six_sections_populated_for_full_template(self) -> None:
        hedz = self.templates["synthetic_hedz"]
        self.assertEqual(len(hedz.sections), 6)
        for i, sec in enumerate(hedz.sections, start=1):
            self.assertEqual(sec["index"], i)
            # Every section in synthetic_hedz is non-empty by construction.
            self.assertNotEqual(sec["type"], 0)

    def test_missing_sections_zero_fill(self) -> None:
        """synthetic_acwow only sets sections 1 and 3 in the fixture;
        sections 2, 4, 5, 6 must default to all-zero entries."""
        acwow = self.templates["synthetic_acwow"]
        self.assertEqual(len(acwow.sections), 6)
        for missing_idx in (2, 4, 5, 6):
            sec = acwow.sections[missing_idx - 1]
            self.assertEqual(sec["index"], missing_idx)
            self.assertEqual(sec["type"], 0)
            self.assertEqual(sec["low_freq"], 0)
            self.assertEqual(sec["low_gain"], 0)
            self.assertEqual(sec["high_freq"], 0)
            self.assertEqual(sec["high_gain"], 0)

    def test_template_first_seen_order_is_stable(self) -> None:
        names = list(self.templates.keys())
        self.assertEqual(names, ["synthetic_hedz", "synthetic_acwow"])

    def test_unrelated_lines_are_ignored(self) -> None:
        """Prose, headings, and bullets that don't match the
        `<name>/filter/...` shape must be silently dropped."""
        noisy = (
            "# Some heading\n"
            "Some prose paragraph.\n"
            "- not_a_filter_bullet: 42\n"
            "- foo/notfilter/something: 5\n"
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
        self.templates = parse_notebooklm_export(
            FIXTURE.read_text(encoding="utf-8")
        )
        self.inventory = to_inventory(self.templates, extracted_from="fixture")

    def test_inventory_format_and_count(self) -> None:
        self.assertEqual(self.inventory["format"], INVENTORY_FORMAT)
        self.assertEqual(self.inventory["template_count"], 2)
        self.assertEqual(self.inventory["extracted_from"], "fixture")

    def test_inventory_round_trips_through_json(self) -> None:
        blob = json.dumps(self.inventory)
        decoded = json.loads(blob)
        self.assertEqual(decoded["format"], INVENTORY_FORMAT)
        hedz = next(t for t in decoded["templates"] if t["name"] == "synthetic_hedz")
        self.assertEqual(hedz["sections"][0]["type"], 3)
        self.assertEqual(hedz["sections"][0]["high_freq"], 116)
        self.assertEqual(hedz["sections"][0]["high_gain"], 87)

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
            self.assertEqual(len(templates), 2)

            names = {t.name for t in templates}
            self.assertEqual(names, {"synthetic_hedz", "synthetic_acwow"})

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
        templates = extract_directory(FIXTURE)
        self.assertIn("synthetic_hedz", templates)
        self.assertIn("synthetic_acwow", templates)

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
