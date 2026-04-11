"""Flagship body generator — E-mu P2K pipeline, 4 shipping bodies.

6 active stages per body, authored in the heritage (type, freq, gain) 0-127
vocabulary. 4 independent corners compiled through the firmware recipes.
Gates filter results. Survivors go to vault.

The pipeline:
  1. Roll 6 stages × 4 corners from body-specific integer ranges
  2. Compile through heritage_coeffs (type 1/2/3 firmware recipes)
  3. 6 active + 6 passthrough per corner (hardware topology)
  4. Gate: peak, dynamic range, talkingness, morph movement
  5. Export compiled-v1 JSON

Usage:
    python tools/flagship_gen.py --body small_talk --count 500
    python tools/flagship_gen.py --body all --count 300
"""
from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

import numpy as np

from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.corner import CornerName
from pyruntime.designer_compile import (
    DesignerCell, make_four_corner_template, compile_four_corner_to_body,
)
from pyruntime.freq_response import cascade_response_db, freq_points

FREQS = freq_points()
_VOCAL_MASK = (FREQS >= 200) & (FREQS <= 5000)
VAULT_DIR = Path(__file__).parent.parent / "vault"


# =============================================================================
# Stage range — one authored stage across 4 corners
# =============================================================================

class SR6:
    """One of 6 stages. Each corner gets independent (type, freq, gain) ranges."""

    def __init__(self, a, b=None, c=None, d=None):
        """Each arg is (type, freq_lo, freq_hi, gain_lo, gain_hi).
        b defaults to a, c defaults to a, d defaults to c."""
        self.a = a
        self.b = b or a
        self.c = c or a
        self.d = d or self.c

    def roll(self, rng: random.Random) -> dict[str, DesignerCell]:
        def _pick(spec):
            t, fl, fh, gl, gh = spec
            return DesignerCell(type=t, freq=rng.randint(fl, fh), gain=rng.randint(gl, gh))
        return {"M0_Q0": _pick(self.a), "M0_Q100": _pick(self.b),
                "M100_Q0": _pick(self.c), "M100_Q100": _pick(self.d)}


class BodySeed:
    def __init__(self, name, key, stages,
                 gate_talk_mid=0.0, gate_talk_min=0.0,
                 gate_peak=48.0, gate_ident=0.20):
        self.name = name
        self.key = key
        self.stages = stages  # exactly 6 SR6
        self.gate_talk_mid = gate_talk_mid
        self.gate_talk_min = gate_talk_min
        self.gate_peak = gate_peak
        self.gate_ident = gate_ident


# =============================================================================
# freq reference (type 1, gain=64):
#   0→72  10→105  20→152  30→220  40→320  50→463  60→686
#   70→993  80→1440  90→2095  100→3066  105→3800  110→4632
#   115→5600  120→7045  125→9200  127→10259
#
# type 1: bandpass resonator      type 2: resonant formant
# type 3: shaped asymmetric       type 0: bypass
# gain 64 = unity. <64 atten. >64 boost.
# =============================================================================

# --- SPEAKER KNOCKERZ ---
# Sub anchor (locked) + character that wakes from dormant to brutal.
# Sub never moves. Character stages: low gain at rest, high gain+freq at morph.
# Q compresses character toward low mids.
SPEAKER_KNOCKERZ = BodySeed(
    "Speaker Knockerz", "speaker_knockerz",
    [
        # S0: sub anchor — locked across all corners
        #          (type, freq_lo, freq_hi, gain_lo, gain_hi)
        SR6(       (1,  0,  5,  66, 78),  # A: ~72-80 Hz, hot
            b=     (1,  0,  5,  68, 80),  # B: Q boosts sub
            c=     (1,  0,  5,  66, 78),  # C: same
            d=     (1,  0,  5,  68, 80)), # D: same
        # S1: sub harmonic
        SR6(       (1,  5, 15,  62, 74),
            b=     (1,  5, 15,  64, 76),
            c=     (1,  5, 15,  62, 74),
            d=     (1,  5, 15,  64, 76)),
        # S2: choke/rip — dormant at rest, wakes with morph
        SR6(       (1, 15, 30,  30, 48),  # A: 130-220 Hz, quiet
            b=     (1, 15, 30,  32, 50),  # B: slightly louder with Q
            c=     (1, 35, 55,  64, 82),  # C: 260-500 Hz, loud — the rip
            d=     (2, 30, 50,  66, 84)), # D: formant type at Q
        # S3: rattle/cry
        SR6(       (1, 25, 45,  28, 44),  # A: 180-360 Hz, quiet
            b=     (1, 25, 45,  30, 46),
            c=     (2, 55, 75,  62, 80),  # C: 500-1000 Hz, the cry
            d=     (2, 50, 70,  64, 82)),
        # S4: stress zone
        SR6(       (2, 35, 55,  24, 40),  # A: 260-500 Hz, barely there
            b=     (2, 35, 55,  26, 42),
            c=     (1, 70, 95,  60, 78),  # C: 993-2400 Hz, stress
            d=     (1, 65, 90,  62, 80)),
        # S5: fracture/air
        SR6(       (1, 50, 70,  20, 36),  # A: 463-993 Hz, ghost
            b=     (1, 50, 70,  22, 38),
            c=     (3, 85,110,  58, 76),  # C: 1700-4632 Hz, shaped fracture
            d=     (3, 80,105,  60, 78)),
    ],
    gate_peak=48.0, gate_ident=0.15,
)

# --- ALUMINUM SIDING ---
# Permanent midrange void. All character in extreme HF.
# Q compresses HF toward center.
ALUMINUM_SIDING = BodySeed(
    "Aluminum Siding", "aluminum_siding",
    [
        # S0: void — midrange scoop, always present
        SR6(       (1, 60, 80,  18, 34),  # 686-1440 Hz, suppressed
            c=     (1, 60, 80,  16, 32)),
        # S1: low warmth floor
        SR6(       (2, 15, 35,  56, 68),  # 130-260 Hz, warm formant
            c=     (2, 15, 35,  54, 66)),
        # S2-S5: HF character — shimmer at rest, tear with morph
        SR6(       (1,100,112,  52, 66),  # A: 3066-5000 Hz
            b=     (1, 95,108,  54, 68),  # B: Q pulls down
            c=     (1,108,118,  62, 78),  # C: 4200-6400 Hz, sharper
            d=     (1,104,115,  64, 80)),
        SR6(       (1,105,116,  48, 62),
            b=     (1,100,112,  50, 64),
            c=     (3,112,122,  60, 76),  # shaped metallic
            d=     (3,108,118,  62, 78)),
        SR6(       (1,110,120,  44, 58),
            b=     (1,106,116,  46, 60),
            c=     (1,118,125,  58, 74),
            d=     (1,114,122,  60, 76)),
        SR6(       (1,115,125,  40, 54),
            b=     (1,110,120,  42, 56),
            c=     (1,122,127,  56, 72),  # near-Nyquist tear
            d=     (1,118,127,  58, 74)),
    ],
    gate_peak=48.0, gate_ident=0.10,
)

# --- SMALL TALK Ah-Ee ---
# Vocal formant morph. F1 drops, F2 rises. Type 2 dominant.
# Q compresses toward nasal center (~1200 Hz).
SMALL_TALK = BodySeed(
    "Small Talk Ah-Ee", "small_talk",
    [
        # S0: F1 — open "Ah" → tight "Ee"
        SR6(       (2, 55, 68,  78, 96),  # A: 500-850 Hz, hot formant
            b=     (2, 50, 62,  80, 98),  # B: Q compresses toward 1200
            c=     (2, 30, 45,  78, 96),  # C: 220-360 Hz — F1 drops
            d=     (2, 35, 50,  80, 98)),
        # S1: F2 — low "Ah" → high "Ee"
        SR6(       (2, 75, 85,  78, 96),  # A: 1100-1700 Hz
            b=     (2, 72, 82,  80, 98),
            c=     (2, 90,102,  78, 96),  # C: 2095-3200 Hz — F2 rises
            d=     (2, 85, 98,  80, 98)),
        # S2: F3 — brightness, slight rise
        SR6(       (2, 95,105,  72, 88),  # A: 2400-3800 Hz
            b=     (2, 90,100,  74, 90),
            c=     (2,100,110,  72, 88),  # C: 3066-4632 Hz
            d=     (2, 95,105,  74, 90)),
        # S3: F4/nasality
        SR6(       (1,105,115,  68, 84),  # A: 3800-5600 Hz
            b=     (1,100,110,  70, 86),
            c=     (1,108,118,  68, 84),
            d=     (1,103,113,  70, 86)),
        # S4: warmth floor — locked
        SR6(       (1,  0, 12,  68, 82),  # 72-115 Hz
            c=     (1,  0, 12,  68, 82)),
        # S5: throat color
        SR6(       (2, 40, 58,  66, 82),  # A: 320-570 Hz
            b=     (2, 45, 62,  68, 84),
            c=     (2, 60, 78,  66, 82),  # C: 686-1200 Hz
            d=     (2, 55, 72,  68, 84)),
    ],
    gate_talk_mid=0.0, gate_talk_min=0.0,
    gate_peak=48.0, gate_ident=0.20,
)

# --- CUL-DE-SAC ---
# Pipe cluster at rest → scattered comb at morph.
# Hum anchor locked. Q compresses the scatter.
CUL_DE_SAC = BodySeed(
    "Cul-De-Sac", "cul_de_sac",
    [
        # S0: hum anchor — locked, always present
        SR6(       (1,  0,  8,  68, 84),  # ~72-100 Hz, strong
            b=     (1,  0,  8,  70, 86),
            c=     (1,  0,  8,  68, 84),
            d=     (1,  0,  8,  70, 86)),
        # S1-S5: pipe at rest, comb at morph
        # At M0: tight cluster 220-600 Hz (type 1/2)
        # At M100: scattered across spectrum (different freqs per stage)
        SR6(       (1, 28, 42,  62, 78),  # A: 200-310 Hz — pipe
            b=     (2, 28, 42,  64, 80),
            c=     (1, 72, 88,  60, 76),  # C: 1000-1600 Hz — shard 1
            d=     (1, 68, 84,  62, 78)),
        SR6(       (2, 32, 46,  62, 78),  # A: 230-370 Hz
            b=     (2, 32, 46,  64, 80),
            c=     (1, 88,102,  58, 74),  # C: 1600-3200 Hz — shard 2
            d=     (1, 84, 98,  60, 76)),
        SR6(       (1, 35, 50,  60, 76),  # A: 260-463 Hz
            b=     (2, 35, 50,  62, 78),
            c=     (2,100,112,  56, 72),  # C: 3066-5000 Hz — shard 3
            d=     (2, 96,108,  58, 74)),
        SR6(       (1, 38, 54,  58, 74),  # A: 280-500 Hz
            b=     (1, 38, 54,  60, 76),
            c=     (1,108,120,  54, 70),  # C: 4200-7045 Hz — shard 4
            d=     (1,104,116,  56, 72)),
        SR6(       (2, 42, 58,  56, 72),  # A: 310-570 Hz
            b=     (2, 42, 58,  58, 74),
            c=     (3,116,127,  52, 68),  # C: 5600-10259 Hz — shard 5
            d=     (3,112,124,  54, 70)),
    ],
    gate_peak=48.0, gate_ident=0.10,
)

ALL_SEEDS = {
    "speaker_knockerz": SPEAKER_KNOCKERZ,
    "aluminum_siding": ALUMINUM_SIDING,
    "small_talk": SMALL_TALK,
    "cul_de_sac": CUL_DE_SAC,
}


# =============================================================================
# Generator + Gates
# =============================================================================

def generate_body(seed: BodySeed, rng: random.Random) -> Body:
    corners = {"M0_Q0": [], "M0_Q100": [], "M100_Q0": [], "M100_Q100": []}
    for sr in seed.stages:
        cells = sr.roll(rng)
        for k in corners:
            corners[k].append(cells[k])
    template = make_four_corner_template(seed.name, corners)
    return compile_four_corner_to_body(template)


def _talk(db):
    vdb = db[_VOCAL_MASK]
    if len(vdb) < 3: return 0.0
    m = float(np.mean(db))
    p = sum(1 for i in range(1, len(vdb)-1)
            if vdb[i] > vdb[i-1] and vdb[i] > vdb[i+1] and vdb[i] > m + 3)
    return min(1.0, p / 4.0)


def _ident(a, b):
    an, bn = a - np.mean(a), b - np.mean(b)
    d = float(np.sqrt(np.sum(an**2) * np.sum(bn**2)))
    return max(0.0, float(np.sum(an * bn) / d)) if d > 1e-10 else 0.0


def _cent(db):
    lin = 10.0 ** (db / 20.0)
    d = float(np.sum(lin))
    return float(np.sum(FREQS * lin) / d) if d > 0 else 1000.0


def gate(body, seed):
    dbs = [cascade_response_db(body.corners.interpolate(m, 0.5), FREQS, SR)
           for m in (0, .25, .5, .75, 1)]

    for i, db in enumerate(dbs):
        if float(np.max(db)) > seed.gate_peak:
            return False, f"peak={np.max(db):.0f}"

    if float(np.max(dbs[2]) - np.min(dbs[2])) < 5:
        return False, "dead"

    talks = [_talk(db) for db in dbs]
    if seed.gate_talk_mid > 0 and talks[2] < seed.gate_talk_mid:
        return False, f"talk_mid={talks[2]:.2f}"
    if seed.gate_talk_min > 0 and min(talks) < seed.gate_talk_min:
        return False, f"talk_min={min(talks):.2f}"

    if _ident(dbs[0], dbs[-1]) < seed.gate_ident:
        return False, "identity"

    if float(np.sqrt(np.mean((dbs[0] - dbs[-1])**2))) < 1.5:
        return False, "static"

    return True, "ok"


def desc(body):
    dbs = [cascade_response_db(body.corners.interpolate(m, 0.5), FREQS, SR) for m in (0, .5, 1)]
    t = [_talk(d) for d in dbs]
    c = [_cent(d) for d in dbs]
    return (f"talk={t[0]:.2f}/{t[1]:.2f}/{t[2]:.2f}  "
            f"cent={c[0]:.0f}→{c[2]:.0f}  peak={max(np.max(d) for d in dbs):.1f}")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--body", default="all")
    parser.add_argument("--count", type=int, default=300)
    parser.add_argument("--out", type=str, default=None)
    args = parser.parse_args()

    keys = list(ALL_SEEDS.keys()) if args.body == "all" else [args.body]

    for key in keys:
        seed = ALL_SEEDS[key]
        out = Path(args.out) if args.out else VAULT_DIR / "_flagship_gen" / key
        out.mkdir(parents=True, exist_ok=True)

        print(f"\n{'='*60}\n  {seed.name} — {args.count} candidates\n{'='*60}")

        ok_n, fails = 0, {}
        for i in range(args.count):
            body = generate_body(seed, random.Random(i * 7919 + hash(key) % 10000))
            passed, reason = gate(body, seed)
            if not passed:
                fails[reason.split("=")[0]] = fails.get(reason.split("=")[0], 0) + 1
                continue

            ok_n += 1
            name = f"{key}_{ok_n:03d}"
            body.name = name
            with open(out / f"{name}.json", "w") as f:
                f.write(body.to_compiled_json(provenance="flagship-gen"))
            print(f"  {name:<30} {desc(body)}")

        print(f"\n  {ok_n}/{args.count} passed")
        for r, c in sorted(fails.items(), key=lambda x: -x[1]):
            print(f"    {r}: {c}")


if __name__ == "__main__":
    main()
