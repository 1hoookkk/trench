"""Author Bismuth Shrapnel — optimizer-driven zero placement.

Different mental: don't hand-author zeros from sonic tables.
Define what Bismuth Shrapnel SOUNDS LIKE as a fitness function,
then search the zero-parameter space for the best configuration.

The fitness function rewards:
  - Dense comb teeth in 3-15 kHz (crystalline HF)
  - Deep void at 800-1200 Hz (permanent scoop)
  - Midpoint audit PASS (stability)
  - More morph motion (expressiveness)

Heritage stages 0-5: P2k_026 poles untouched.
Character stages 6-11: poles fixed at interleaved HF positions,
ZEROS found by random search over 1000 candidates.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]

import sys
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.corner import CornerName, CornerState, CornerArray
from pyruntime.stage_math import resonator_with_zero
from pyruntime.stage_params import StageParams
from pyruntime.constants import SR
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.analysis import midpoint_audit

SOURCE = ROOT / "cartridges" / "p2k" / "P2k_026.json"
OUT_DIR = ROOT / "cartridges" / "candidates"

FREQS = freq_points(n=512, sr=SR)

# Character poles: interleaved with Millennium heritage
COMB_POLES = [2700, 4900, 7300, 9400, 11400, 13500]

# Search bounds for zeros (Hz) — anywhere in 200-18000 Hz
ZERO_LO = 200.0
ZERO_HI = 18000.0
# Search bounds for zero radius — 0.05 to 0.70
ZR_LO = 0.05
ZR_HI = 0.70

# Fixed character stage parameters
CHAR_VAL1 = 0.0   # unity gain

# Corner r values (fixed — only zeros are searched)
CORNER_R = {
    CornerName.A: 0.30,   # M0_Q0
    CornerName.B: 0.35,   # M0_Q100
    CornerName.C: 0.45,   # M100_Q0
    CornerName.D: 0.12,   # M100_Q100 (near passthrough)
}


def build_candidate(base: Body, zero_config: np.ndarray) -> Body:
    """Build a body from a zero configuration vector.

    zero_config: shape (6, 2) — [zero_hz, zero_r] per character stage.
    Same zeros at all active corners (only r scales by corner).
    """
    new_corners = {}
    for cn in CornerName:
        src = base.corners.corner(cn)
        char_r = CORNER_R[cn]

        stages = list(src.stages[:6])

        for i in range(6):
            zero_hz = float(zero_config[i, 0])
            zero_r_base = float(zero_config[i, 1])
            # Scale zero_r by corner r ratio
            zero_r = zero_r_base * (char_r / 0.45)
            zero_r = max(0.01, min(0.70, zero_r))
            sp = resonator_with_zero(COMB_POLES[i], char_r, CHAR_VAL1, zero_hz, zero_r)
            stages.append(sp)

        while len(stages) < 12:
            stages.append(StageParams.passthrough())

        new_corners[cn] = CornerState(stages=stages, boost=src.boost)

    return Body(
        name="Bismuth Shrapnel candidate",
        corners=CornerArray(
            a=new_corners[CornerName.A],
            b=new_corners[CornerName.B],
            c=new_corners[CornerName.C],
            d=new_corners[CornerName.D],
        ),
        boost=base.boost,
    )


def _count_peaks(db: np.ndarray, freqs: np.ndarray, lo: float, hi: float, threshold: float = 3.0) -> int:
    """Count peaks above mean+threshold in a frequency band."""
    mask = (freqs >= lo) & (freqs <= hi)
    band = db[mask]
    if len(band) < 3:
        return 0
    mean = float(np.mean(band))
    count = 0
    for i in range(1, len(band) - 1):
        if band[i] > band[i-1] and band[i] > band[i+1] and band[i] > mean + threshold:
            count += 1
    return count


def bismuth_fitness(body: Body) -> float:
    """Score a body on Bismuth Shrapnel criteria. Higher = better.

    Rewards:
      - Dense HF comb teeth at M0_Q0 (the Bismuth identity)
      - Deep 1 kHz void (permanent scoop)
      - Positive spectral tilt (HF > LF)
      - Open state NOT silent (M100 must still have energy)
      - Morph motion
    Penalties:
      - Midpoint audit failure (hard kill)
      - Open state peak below -10 dB (too quiet)
      - Void above 0 dB (not scooped enough)
    """
    score = 0.0

    # Hard gate: midpoint audit
    audit = midpoint_audit(body, peak_limit_db=35.0)
    if not audit["passed"]:
        return -1000.0 + (35.0 - audit["worst_peak_db"])

    # --- M0_Q0: Dull Silver (the comb state) ---
    enc_m0 = body.corners.interpolate(0.0, 0.0)
    db_m0 = cascade_response_db(enc_m0, FREQS, SR)

    # Comb teeth in 3-15 kHz
    teeth = _count_peaks(db_m0, FREQS, 3000, 15000, threshold=3.0)
    score += teeth * 8.0

    # HF energy (mean dB in 3-15 kHz) — reward positive
    mask_hf = (FREQS >= 3000) & (FREQS <= 15000)
    hf_energy = float(np.mean(db_m0[mask_hf]))
    score += hf_energy * 1.5

    # 1 kHz void depth (mean dB in 800-1200 Hz)
    mask_void = (FREQS >= 800) & (FREQS <= 1200)
    void_db = float(np.mean(db_m0[mask_void]))
    void_depth = hf_energy - void_db  # positive = void is below HF
    score += max(0, void_depth) * 2.0

    # Spectral tilt: HF minus LF — reward positive tilt
    mask_lf = (FREQS >= 20) & (FREQS <= 500)
    lf_energy = float(np.mean(db_m0[mask_lf]))
    tilt = hf_energy - lf_energy
    score += max(0, tilt) * 1.0

    # --- M100_Q0: open state must not be dead ---
    enc_m1 = body.corners.interpolate(1.0, 0.0)
    db_m1 = cascade_response_db(enc_m1, FREQS, SR)
    open_peak = float(np.max(db_m1))
    if open_peak < -10.0:
        score -= (abs(open_peak) - 10) * 3.0  # heavy penalty for dead open state
    else:
        score += open_peak * 0.5  # mild reward for audible open state

    # --- M100_Q100: should still have some life ---
    enc_m1q1 = body.corners.interpolate(1.0, 1.0)
    db_m1q1 = cascade_response_db(enc_m1q1, FREQS, SR)
    m1q1_peak = float(np.max(db_m1q1))
    if m1q1_peak < -5.0:
        score -= (abs(m1q1_peak) - 5) * 2.0

    # --- Morph distance ---
    morph_dist = float(np.sqrt(np.mean((db_m0 - db_m1) ** 2)))
    score += morph_dist * 0.5

    # --- M0.5 midpoint should have HF character too ---
    enc_mid = body.corners.interpolate(0.5, 0.0)
    db_mid = cascade_response_db(enc_mid, FREQS, SR)
    mid_teeth = _count_peaks(db_mid, FREQS, 3000, 15000, threshold=2.0)
    score += mid_teeth * 3.0

    return score


def random_search(base: Body, n_candidates: int = 500) -> tuple[np.ndarray, float]:
    """Search the zero-parameter space for the best Bismuth configuration."""
    rng = np.random.default_rng(42)
    best_config = None
    best_score = -9999.0

    for i in range(n_candidates):
        # Random zero config: 6 stages, each (zero_hz, zero_r)
        zero_hz = rng.uniform(ZERO_LO, ZERO_HI, size=6)
        zero_r = rng.uniform(ZR_LO, ZR_HI, size=6)
        config = np.stack([zero_hz, zero_r], axis=1)

        body = build_candidate(base, config)
        score = bismuth_fitness(body)

        if score > best_score:
            best_score = score
            best_config = config.copy()
            if (i + 1) % 100 == 0 or i == 0:
                print(f"  [{i+1}/{n_candidates}] best={best_score:.1f} "
                      f"teeth/void/morph components")

    return best_config, best_score


def main() -> None:
    base = Body.from_json(str(SOURCE))

    print(f"Searching 2000 random zero configurations...")
    t0 = time.time()
    best_config, best_score = random_search(base, n_candidates=2000)
    elapsed = time.time() - t0
    print(f"\nSearch complete in {elapsed:.1f}s")
    print(f"Best score: {best_score:.1f}")
    print(f"\nBest zero config (per character stage):")
    for i in range(6):
        print(f"  S{6+i} pole={COMB_POLES[i]}Hz → zero={best_config[i,0]:.0f}Hz "
              f"zero_r={best_config[i,1]:.3f}")

    # Build and save the winner
    body = build_candidate(base, best_config)
    body.name = "Bismuth Shrapnel v1"

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    kf_path = OUT_DIR / "Bismuth_Shrapnel_v1.keyframe.json"
    cc_path = OUT_DIR / "Bismuth_Shrapnel_v1.compiled.json"

    raw = body.to_json()
    d = json.loads(raw)
    d["stages"] = 12
    kf_path.write_text(json.dumps(d, indent=2), encoding="utf-8")
    cc_path.write_text(
        body.to_compiled_json(provenance="optimizer-search-v1"),
        encoding="utf-8",
    )
    print(f"\nwrote {kf_path}")
    print(f"wrote {cc_path}")


if __name__ == "__main__":
    main()
