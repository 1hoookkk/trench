#!/usr/bin/env python3
"""
build_filter_reference.py — Build machine-readable E-MU Z-plane filter reference.

Reads:
  - resources/Filter/*.xml              (71 MorphDesigner user presets)
  - cartridges/engine/_source/heritage_designer_sections.json (83 factory P2K templates)

Outputs:
  - resources/emu_zplane_reference.xml  (enriched reference for LLM/NotebookLM)

Includes raw section params + compiled frequency/radius/kernel at morph endpoints.
Stdlib only — no numpy, no pyruntime imports.
"""

import json
import math
import struct
import sys
import xml.etree.ElementTree as ET
from pathlib import Path
from xml.dom import minidom

REPO = Path(__file__).resolve().parent.parent
FILTER_DIR = REPO / "resources" / "Filter"
HERITAGE_JSON = REPO / "cartridges" / "engine" / "_source" / "heritage_designer_sections.json"
OUTPUT = REPO / "resources" / "emu_zplane_reference.xml"

# ── Firmware compile math (from bake_hedz_const.py) ──────────────────────

SR = 44100
SR_FAMILY = {44100: 0, 48000: 1, 96000: 2, 192000: 3}
FW_BASE = [18, 18, 4, 1]
FW_SCALE = [220, 220, 200, 177]
COMBINE_K = 4.0

SECTION_TYPE_NAMES = {0: "off", 1: "lowpass/highpass", 2: "EQ (parametric)", 3: "lowpass (alt)"}


def f32(x):
    return struct.unpack("<f", struct.pack("<f", x))[0]


def mf_decode(packed):
    packed &= 0xFFFF
    if packed == 0xFFFF:
        return 1.0
    if packed == 0x0000:
        return 0.0
    u = packed + 1
    exp = ((u >> 12) & 0xF) - 15
    mant = u & 0xFFF
    if exp < -14:
        return mant / 134217728.0
    scaled = (mant | 0x1000) / 8192.0
    return scaled * (2.0 ** exp)


def fw_freq_value(packed_freq, sr=SR):
    idx = SR_FAMILY.get(sr, 0)
    return (FW_SCALE[idx] * packed_freq) // 128 + FW_BASE[idx]


def fw_radius(freq_value):
    return (freq_value * 124) // 256 + 118


def fw_gain_offset(packed_gain, shift=0):
    raw = (packed_gain - 64) // 2 + shift
    return max(-32, min(31, raw))


def clamp255(v):
    return max(0, min(255, v))


def fw_words_to_kernel(w0, w1, w2, w3, w4):
    d = [mf_decode(w) for w in [w0, w1, w2, w3, w4]]
    return (
        f32(d[0] * COMBINE_K + d[1]),
        f32(d[1]),
        f32(d[2] * COMBINE_K + d[3]),
        f32(d[3]),
        f32(d[4]),
    )


def type1_kernel(fp, gp, shift, sr):
    fv = fw_freq_value(fp, sr)
    rad = fw_radius(fv)
    go = fw_gain_offset(gp, shift)
    w0 = fv << 8
    w1 = clamp255(rad + go) << 8
    w2 = fv << 8
    w3 = clamp255(rad - go) << 8
    w4 = 0xE000
    return fw_words_to_kernel(w0, w1, w2, w3, w4)


def type2_kernel(fp, gp, shift, sr):
    idx = SR_FAMILY.get(sr, 0)
    fv = fw_freq_value(fp, sr)
    rad = fw_radius(fv)
    go = fw_gain_offset(gp, shift)
    if idx < 2:
        w0 = 0xEC << 8
        w1 = 0xFF << 8
        w4_val = fv + 0xF5
    else:
        w0 = 0xE1 << 8
        w1 = 0xF0 << 8
        w4_val = fv
    w2 = fv << 8
    w3 = clamp255(rad - go) << 8
    w4 = w4_val << 8
    return fw_words_to_kernel(w0, w1, w2, w3, w4)


def type3_freq_compression(freq_value, shift):
    if freq_value > 0xDB and shift < 0:
        freq_value = (((freq_value - 220) * (shift + 32)) >> 5) + 220
    return freq_value


def type3_kernel(fp, gp, shift, sr):
    idx = SR_FAMILY.get(sr, 0)
    fv = fw_freq_value(fp, sr)
    go = fw_gain_offset(gp, shift)
    fvc = type3_freq_compression(fv, shift)
    rad = fw_radius(fvc)
    base = FW_BASE[idx]
    w0 = base << 8
    w1 = ((base * 124) // 256 + 150) << 8
    w2 = fvc << 8
    w3 = clamp255(rad - go) << 8
    if idx < 2:
        c4_raw = (fvc - 18) * (-12) + (-8192)
        w4 = c4_raw & 0xFFFF
    else:
        w4 = 0xE000
    return fw_words_to_kernel(w0, w1, w2, w3, w4)


def compile_stage(type_id, freq_packed, gain_packed, shift=0, sr=SR):
    if type_id <= 0:
        return (1.0, 0.0, 0.0, 0.0, 0.0)
    if type_id == 1:
        return type1_kernel(freq_packed, gain_packed, shift, sr)
    if type_id == 2:
        return type2_kernel(freq_packed, gain_packed, shift, sr)
    if type_id == 3:
        return type3_kernel(freq_packed, gain_packed, shift, sr)
    return type1_kernel(freq_packed, gain_packed, shift, sr)


def pair_freq_hz(a1, a2, sr=SR):
    """Frequency from a conjugate pair encoded as DF2T coefficients.
    Given denominator 1 + a1*z^-1 + a2*z^-2 (or numerator c0 + c1*z^-1 + c2*z^-2):
      r = sqrt(a2), theta = acos(-a1 / (2*r)), freq = theta * sr / (2*pi).
    Returns (freq_hz, radius)."""
    if a2 <= 0:
        return (0.0, 0.0)
    r = math.sqrt(a2)
    if r < 1e-12:
        return (0.0, 0.0)
    cos_theta = -a1 / (2.0 * r)
    cos_theta = max(-1.0, min(1.0, cos_theta))
    theta = math.acos(cos_theta)
    return (theta * sr / (2.0 * math.pi), r)


# ── Section analysis ─────────────────────────────────────────────────────

def analyze_section(type_id, freq, gain, shift=0, sr=SR):
    """Return a dict describing one compiled section."""
    fv = fw_freq_value(freq, sr)
    rad = fw_radius(fv)
    kernel = compile_stage(type_id, freq, gain, shift, sr)
    c0, c1, c2, c3, c4 = kernel

    # Extract zero (numerator) and pole (denominator) frequencies.
    # H(z) = (c0 + c1*z^-1 + c2*z^-2) / (1 + c3*z^-1 + c4*z^-2)
    zero_hz, zero_r = pair_freq_hz(c1, c2, sr)
    pole_hz, pole_r = pair_freq_hz(c3, c4, sr)

    return {
        "type_id": type_id,
        "type_name": SECTION_TYPE_NAMES.get(type_id, f"unknown({type_id})"),
        "freq_packed": freq,
        "gain_packed": gain,
        "freq_value": fv,
        "radius": rad,
        "zero_hz": round(zero_hz, 2),
        "zero_radius": round(zero_r, 6),
        "pole_hz": round(pole_hz, 2),
        "pole_radius": round(pole_r, 6),
        "kernel": [round(c, 8) for c in kernel],
    }


# ── Parse XML templates ─────────────────────────────────────────────────

def parse_xml_template(path):
    """Parse a resources/Filter/*.xml file into a dict."""
    tree = ET.parse(path)
    root = tree.getroot()
    name = root.get("name", path.stem)

    filt = root.find("filter")
    frequency = float(filt.findtext("frequency", "0"))
    gain = float(filt.findtext("gain", "0"))
    type_absolute = int(filt.findtext("type-absolute", "108"))

    sections = []
    for sec in filt.findall("designer-section"):
        sections.append({
            "index": int(sec.get("index", 0)),
            "type": int(sec.findtext("type", "0")),
            "low_freq": int(sec.findtext("low-freq", "0")),
            "low_gain": int(sec.findtext("low-gain", "0")),
            "high_freq": int(sec.findtext("high-freq", "0")),
            "high_gain": int(sec.findtext("high-gain", "0")),
        })

    return {
        "name": name,
        "source": str(path.relative_to(REPO)),
        "origin": "user-preset",
        "frequency": frequency,
        "gain": gain,
        "type_absolute": type_absolute,
        "sections": sorted(sections, key=lambda s: s["index"]),
    }


def parse_heritage_templates(path):
    """Parse heritage_designer_sections.json into a list of dicts."""
    with open(path) as f:
        data = json.load(f)

    templates = []
    for tmpl in data.get("templates", []):
        sections = []
        for sec in tmpl.get("sections", []):
            sections.append({
                "index": sec["index"],
                "type": sec["type"],
                "low_freq": sec["low_freq"],
                "low_gain": sec["low_gain"],
                "high_freq": sec["high_freq"],
                "high_gain": sec["high_gain"],
            })
        templates.append({
            "name": tmpl["name"],
            "source": str(path.relative_to(REPO)),
            "origin": "factory-p2k",
            "type_absolute": tmpl.get("type_absolute", 0),
            "frequency": tmpl.get("frequency", 0),
            "gain": tmpl.get("gain", 0),
            "shape_id": tmpl.get("shape_id"),
            "sections": sorted(sections, key=lambda s: s["index"]),
        })
    return templates


# ── Build enriched reference ─────────────────────────────────────────────

def build_filter_entry(tmpl):
    """Analyze a template at morph=0, morph=0.5, morph=1."""
    entry = {
        "name": tmpl["name"],
        "source": tmpl["source"],
        "origin": tmpl["origin"],
        "type_absolute": tmpl.get("type_absolute", 0),
        "sections": tmpl["sections"],
        "morph_endpoints": {},
    }

    active_count = sum(1 for s in tmpl["sections"] if s["type"] > 0)
    entry["active_sections"] = active_count
    entry["order"] = active_count * 2

    for morph_label, morph_val in [("morph_0", 0.0), ("morph_50", 0.5), ("morph_100", 1.0)]:
        stages = []
        for sec in tmpl["sections"]:
            freq = round(sec["low_freq"] + morph_val * (sec["high_freq"] - sec["low_freq"]))
            gain = round(sec["low_gain"] + morph_val * (sec["high_gain"] - sec["low_gain"]))
            analysis = analyze_section(sec["type"], freq, gain)
            stages.append(analysis)
        entry["morph_endpoints"][morph_label] = stages

    return entry


def entry_to_xml(entry):
    """Convert one filter entry to XML element."""
    el = ET.Element("filter")
    el.set("name", entry["name"])
    el.set("origin", entry["origin"])
    el.set("active-sections", str(entry["active_sections"]))
    el.set("order", str(entry["order"]))

    src = ET.SubElement(el, "source")
    src.text = entry["source"]

    if entry.get("type_absolute"):
        ta = ET.SubElement(el, "type-absolute")
        ta.text = str(entry["type_absolute"])

    raw = ET.SubElement(el, "raw-sections")
    for sec in entry["sections"]:
        s = ET.SubElement(raw, "section")
        s.set("index", str(sec["index"]))
        s.set("type", str(sec["type"]))
        s.set("type-name", SECTION_TYPE_NAMES.get(sec["type"], "unknown"))
        s.set("low-freq", str(sec["low_freq"]))
        s.set("low-gain", str(sec["low_gain"]))
        s.set("high-freq", str(sec["high_freq"]))
        s.set("high-gain", str(sec["high_gain"]))

    for morph_label, stages in entry["morph_endpoints"].items():
        mp = ET.SubElement(el, "compiled")
        mp.set("morph", morph_label)
        for i, stage in enumerate(stages):
            st = ET.SubElement(mp, "stage")
            st.set("index", str(i + 1))
            st.set("type", stage["type_name"])
            if stage["type_id"] == 0:
                st.set("status", "off")
                continue
            st.set("freq-packed", str(stage["freq_packed"]))
            st.set("gain-packed", str(stage["gain_packed"]))
            st.set("freq-value", str(stage["freq_value"]))
            st.set("radius", str(stage["radius"]))
            st.set("zero-hz", str(stage["zero_hz"]))
            st.set("zero-radius", str(stage["zero_radius"]))
            st.set("pole-hz", str(stage["pole_hz"]))
            st.set("pole-radius", str(stage["pole_radius"]))
            k = ET.SubElement(st, "kernel")
            k.set("c0", str(stage["kernel"][0]))
            k.set("c1", str(stage["kernel"][1]))
            k.set("c2", str(stage["kernel"][2]))
            k.set("c3", str(stage["kernel"][3]))
            k.set("c4", str(stage["kernel"][4]))

    return el


def pretty_xml(root):
    raw = ET.tostring(root, encoding="unicode", xml_declaration=False)
    dom = minidom.parseString(raw)
    lines = dom.toprettyxml(indent="  ", encoding=None).split("\n")
    # Skip minidom's own xml declaration (first line)
    return "\n".join(line for line in lines[1:] if line.strip())


def main():
    # Collect all templates
    all_templates = []

    # XML user presets
    xml_count = 0
    if FILTER_DIR.is_dir():
        for xml_path in sorted(FILTER_DIR.glob("*.xml")):
            tmpl = parse_xml_template(xml_path)
            all_templates.append(tmpl)
            xml_count += 1

    # Heritage factory templates
    heritage_count = 0
    if HERITAGE_JSON.is_file():
        heritage = parse_heritage_templates(HERITAGE_JSON)
        all_templates.extend(heritage)
        heritage_count = len(heritage)

    print(f"Loaded {xml_count} XML user presets, {heritage_count} factory P2K templates")

    # Build enriched entries
    entries = [build_filter_entry(t) for t in all_templates]

    # ── Freq lookup table: packed 0-127 → Hz per type ──────────────────
    freq_table = {}
    for type_id in [1, 2, 3]:
        rows = []
        for pf in range(128):
            kernel = compile_stage(type_id, pf, 64, 0, SR)
            fv = fw_freq_value(pf, SR)
            zh, zr = pair_freq_hz(kernel[1], kernel[2], SR)
            ph, pr = pair_freq_hz(kernel[3], kernel[4], SR)
            rows.append({"packed": pf, "freq_value": fv,
                          "zero_hz": round(zh, 2), "zero_radius": round(zr, 6),
                          "pole_hz": round(ph, 2), "pole_radius": round(pr, 6)})
        freq_table[type_id] = rows

    # ── Gain sweep table: gain 0-127 effect per type at freq=64 ──────
    gain_table = {}
    for type_id in [1, 2, 3]:
        rows = []
        for gp in range(128):
            kernel = compile_stage(type_id, 64, gp, 0, SR)
            zh, zr = pair_freq_hz(kernel[1], kernel[2], SR)
            ph, pr = pair_freq_hz(kernel[3], kernel[4], SR)
            rows.append({"gain_packed": gp,
                          "zero_hz": round(zh, 2), "zero_radius": round(zr, 6),
                          "pole_hz": round(ph, 2), "pole_radius": round(pr, 6)})
        gain_table[type_id] = rows

    # Build XML document
    root = ET.Element("emu-zplane-filter-reference")
    root.set("generated-by", "tools/build_filter_reference.py")
    root.set("sample-rate", str(SR))
    root.set("total-filters", str(len(entries)))
    root.set("xml-user-presets", str(xml_count))
    root.set("factory-p2k-templates", str(heritage_count))

    meta = ET.SubElement(root, "description")
    meta.text = (
        "Machine-readable reference for E-MU Z-plane MorphDesigner filters. "
        "Each filter entry includes raw section parameters (type, freq, gain per morph endpoint) "
        "and compiled biquad coefficients at morph=0%, 50%, and 100%. "
        "Section types: 0=off, 1=lowpass/highpass, 2=EQ(parametric), 3=lowpass(alt). "
        "Freq and gain are packed integers 0-127. Kernel coefficients are DF2T direct form: "
        "y=c0*x+w1, w1=c1*x-c3*y+w2, w2=c2*x-c4*y. "
        "Six sections cascade serially (12 poles max). "
        "Morph interpolates linearly between low (frame A) and high (frame B) parameters."
    )

    arch = ET.SubElement(root, "architecture")
    arch.text = (
        "The Z-plane filter is a 12-stage serial DF2T biquad cascade. "
        "6 sections × 2 poles each = 12 poles maximum. "
        "Each section can be independently set to off (passthrough), "
        "lowpass/highpass (type 1), EQ parametric (type 2), or lowpass alt (type 3). "
        "The Morph parameter interpolates all 6 sections simultaneously between "
        "two programmed frames (A at morph=0, B at morph=100). "
        "The Q/Gain wheel applies a global offset: for EQ stages it shifts gain, "
        "for LP/HP stages it shifts frequency by ±50%."
    )

    # ── Freq lookup tables ──────────────────────────────────────────────
    freq_ref = ET.SubElement(root, "freq-lookup")
    freq_ref.set("description",
        "Packed freq value (0-127) to exact pole frequency in Hz at 44100 Hz sample rate. "
        "One table per section type. Use this to find what packed value gives a target frequency. "
        "Gain held at 64 (neutral) for isolation.")

    for type_id in [1, 2, 3]:
        tbl = ET.SubElement(freq_ref, "type")
        tbl.set("id", str(type_id))
        tbl.set("name", SECTION_TYPE_NAMES[type_id])
        for row in freq_table[type_id]:
            r = ET.SubElement(tbl, "entry")
            r.set("packed", str(row["packed"]))
            r.set("freq-value", str(row["freq_value"]))
            r.set("zero-hz", str(row["zero_hz"]))
            r.set("zero-radius", str(row["zero_radius"]))
            r.set("pole-hz", str(row["pole_hz"]))
            r.set("pole-radius", str(row["pole_radius"]))

    # ── Gain behavior tables ─────────────────────────────────────────
    gain_ref = ET.SubElement(root, "gain-lookup")
    gain_ref.set("description",
        "Gain value (0-127) effect on pole placement per section type. "
        "Freq held at 64 (midrange) for isolation. Shows how gain reshapes the filter.")

    gain_semantics = {
        1: ("For LP/HP stages, gain controls resonance sharpness via radius split between "
            "conjugate pole pairs. gain=0 → narrow/sharp, gain=64 → neutral, gain=127 → wide/dull. "
            "The big Q Wheel is Q Offset (-50% to +50%)."),
        2: ("For EQ stages, gain controls boost/cut depth via radius asymmetry. "
            "gain=0 → maximum cut, gain=64 → flat/neutral, gain=127 → maximum boost (+24 dB). "
            "The big Gain Wheel is Gain Offset (-24 to +24 dB). Two controls sum and clip."),
        3: ("For LP alt stages, gain controls resonance similarly to type 1 but with "
            "frequency compression at high freq values (>0xDB). gain=0 → sharp, gain=127 → wide. "
            "Includes additional high-frequency rolloff via c4 term."),
    }

    for type_id in [1, 2, 3]:
        tbl = ET.SubElement(gain_ref, "type")
        tbl.set("id", str(type_id))
        tbl.set("name", SECTION_TYPE_NAMES[type_id])
        sem = ET.SubElement(tbl, "behavior")
        sem.text = gain_semantics[type_id]
        for row in gain_table[type_id]:
            r = ET.SubElement(tbl, "entry")
            r.set("gain", str(row["gain_packed"]))
            r.set("zero-hz", str(row["zero_hz"]))
            r.set("zero-radius", str(row["zero_radius"]))
            r.set("pole-hz", str(row["pole_hz"]))
            r.set("pole-radius", str(row["pole_radius"]))

    user_group = ET.SubElement(root, "group")
    user_group.set("name", "MorphDesigner User Presets")
    user_group.set("count", str(xml_count))
    user_group.set("source", "resources/Filter/*.xml")

    factory_group = ET.SubElement(root, "group")
    factory_group.set("name", "Factory P2K Templates")
    factory_group.set("count", str(heritage_count))
    factory_group.set("source", "cartridges/engine/_source/heritage_designer_sections.json")

    for entry in entries:
        parent = user_group if entry["origin"] == "user-preset" else factory_group
        parent.append(entry_to_xml(entry))

    # Write
    xml_str = '<?xml version="1.0" encoding="UTF-8"?>\n' + pretty_xml(root)
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT.write_text(xml_str, encoding="utf-8")
    print(f"Wrote {OUTPUT} ({len(xml_str):,} bytes, {len(entries)} filters)")


if __name__ == "__main__":
    main()
