# Shipping Bodies V1 Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Author 4 shipping bodies (Speaker Knockerz, Aluminum Siding, Small Talk Ah-Ee, Cul-De-Sac) using the P2K StageParams pipeline with trajectory fitness scoring.

**Architecture:** Each body is built from direct pole/zero placement via `stage_math.py` (not the heritage integer compiler, which produces attenuated cascades). Formant targets come from the VOWELS table and acoustic references. Each body gets 4 independent corners (M0_Q0, M0_Q100, M100_Q0, M100_Q100) with the Morpheus comb pattern (alternating peak/notch) adapted to P2K-grade radii (0.975-0.999). Bodies are scored with the trajectory/continuity fitness stack and validated via trench-mcp `analyze_state` and `scan_vault`.

**Tech Stack:** Python (pyruntime), stage_math.py, zero_law.py, encode.py, corner.py, body.py, forge_optimize.py (scoring), trench-mcp (validation)

**Key Reference Files:**
- `pyruntime/stage_math.py` — `resonator()`, `resonator_with_zero()`, `zero_forced()`
- `pyruntime/zero_law.py` — ContourKind, resolve_contour()
- `pyruntime/encode.py` — `raw_to_encoded()` (StageParams → EncodedCoeffs)
- `pyruntime/corner.py` — CornerArray, bilinear interpolation
- `pyruntime/body.py` — Body, JSON I/O
- `pyruntime/forge_optimize.py` — `_eval_body()`, `_extract_formants()`, `_trajectory_score()`
- `cartridges/heritage/Ooh_to_Eee_(approx).json` — reference vocal body (composite=0.96)
- `cartridges/p2k/P2k_013.json` — Talking Hedz P2K (composite=0.64)
- `SHIPPING.md` — curatorial brief and morph path descriptions

**P2K-grade parameter ranges (from working heritage bodies):**
- `radius`: 0.975–0.999 (sharp resonances, +6 to +38 dB peaks)
- `val1`: -0.438 (P2K standard weight, produces ~0.66 weight factor)
- Zero placement: INTERIOR_ZERO with zero_r 0.70-0.75, offset 1-2 octaves from pole
- Sample rate: 39062.5 Hz (authoring domain)

---

### Task 1: Build body authoring recipe script

**Files:**
- Create: `pyruntime/recipes/author_shipping.py`
- Reference: `pyruntime/stage_math.py`, `pyruntime/zero_law.py`, `pyruntime/corner.py`, `pyruntime/body.py`

- [ ] **Step 1: Create the recipe scaffold**

```python
"""Shipping body authoring — 4 flagship bodies for TRENCH V1.

Each body is authored as 4 independent corners using direct pole/zero
placement at P2K-grade radii. NOT the heritage integer compiler.

Run: python -m pyruntime.recipes.author_shipping
"""
import json
import math
import os
from pathlib import Path

from pyruntime.body import Body
from pyruntime.constants import SR, TWO_PI, NUM_BODY_STAGES
from pyruntime.corner import CornerArray, CornerState
from pyruntime.encode import raw_to_encoded
from pyruntime.stage_math import resonator, resonator_with_zero, zero_forced
from pyruntime.stage_params import StageParams

VAL1_P2K = -0.438    # P2K-calibrated weight (0.66)
ZERO_R = 0.72        # Interior zero radius (from calibration data)
OUT_DIR = Path(__file__).parent.parent.parent / "cartridges" / "candidates"


def _interior_zero_stage(freq_hz: float, radius: float,
                         zero_offset_octaves: float = 1.5) -> StageParams:
    """Single stage with interior zero, P2K-calibrated."""
    zero_freq = freq_hz * (2.0 ** zero_offset_octaves)
    zero_freq = max(20.0, min(SR * 0.48, zero_freq))
    return resonator_with_zero(freq_hz, radius, VAL1_P2K, zero_freq, ZERO_R)


def _passthrough() -> StageParams:
    return StageParams.passthrough()


def _build_corner(stages: list[StageParams], boost: float = 4.0) -> CornerState:
    """Pad to NUM_BODY_STAGES and create CornerState."""
    while len(stages) < NUM_BODY_STAGES:
        stages.append(_passthrough())
    return CornerState(stages=stages, boost=boost)


def _build_body(name: str, corners: dict) -> Body:
    """Build body from 4 named corners."""
    ca = CornerArray(
        a=corners["M0_Q0"],
        b=corners["M0_Q100"],
        c=corners["M100_Q0"],
        d=corners["M100_Q100"],
    )
    return Body(name=name, corners=ca, boost=4.0)
```

- [ ] **Step 2: Verify scaffold imports**

Run: `python -c "from pyruntime.recipes.author_shipping import _interior_zero_stage, _build_corner; print('OK')"`
Expected: OK

- [ ] **Step 3: Commit**

```bash
git add pyruntime/recipes/author_shipping.py
git commit -m "feat: shipping body authoring scaffold with P2K-grade stage helpers"
```

---

### Task 2: Author Small Talk Ah-Ee (body 3)

Start with the vocal body — it has the clearest acoustic targets and the best reference (Ooh_to_Eee composite=0.96).

**Files:**
- Modify: `pyruntime/recipes/author_shipping.py`
- Output: `cartridges/candidates/SmallTalk_AhEe_v1.json`

**Acoustic targets (from VOWELS table):**
- Ah: F1=768, F2=1189, F3=2555, F4=3508
- Ee: F1=342, F2=2322, F3=3000, F4=3657

**Construction pattern (from hedz/Talking Hedz calibration):**
- 3 active stages per corner with interior zeros
- Stage 0: F1 formant anchor (r=0.990-0.995)
- Stage 1: F2 formant bite (r=0.985-0.995)
- Stage 2: F3/F4 air (r=0.975-0.985)
- Q axis: tighten radii at Q=1 (sharper formants under pressure)

- [ ] **Step 1: Implement author_small_talk()**

```python
def author_small_talk() -> Body:
    """Body 3: Small Talk Ah-Ee — Ah vowel to Ee vowel morph.

    3 active formant stages with interior zeros.
    Q tightens radii (sharper formants under pressure).
    """
    # -- Formant targets --
    AH = {"f1": 768, "f2": 1189, "f3": 2555}
    EE = {"f1": 342, "f2": 2322, "f3": 3000}

    # -- Corner A: M0_Q0 (Ah, relaxed) --
    a = _build_corner([
        _interior_zero_stage(AH["f1"], 0.990, zero_offset_octaves=1.5),
        _interior_zero_stage(AH["f2"], 0.988, zero_offset_octaves=1.3),
        _interior_zero_stage(AH["f3"], 0.978, zero_offset_octaves=1.0),
    ])

    # -- Corner B: M0_Q100 (Ah, pressured — sharper formants) --
    b = _build_corner([
        _interior_zero_stage(AH["f1"], 0.996, zero_offset_octaves=1.5),
        _interior_zero_stage(AH["f2"], 0.994, zero_offset_octaves=1.3),
        _interior_zero_stage(AH["f3"], 0.985, zero_offset_octaves=1.0),
    ])

    # -- Corner C: M100_Q0 (Ee, relaxed) --
    c = _build_corner([
        _interior_zero_stage(EE["f1"], 0.990, zero_offset_octaves=1.5),
        _interior_zero_stage(EE["f2"], 0.988, zero_offset_octaves=1.3),
        _interior_zero_stage(EE["f3"], 0.978, zero_offset_octaves=1.0),
    ])

    # -- Corner D: M100_Q100 (Ee, pressured) --
    d = _build_corner([
        _interior_zero_stage(EE["f1"], 0.996, zero_offset_octaves=1.5),
        _interior_zero_stage(EE["f2"], 0.994, zero_offset_octaves=1.3),
        _interior_zero_stage(EE["f3"], 0.985, zero_offset_octaves=1.0),
    ])

    return _build_body("SmallTalk_AhEe", {
        "M0_Q0": a, "M0_Q100": b, "M100_Q0": c, "M100_Q100": d,
    })
```

- [ ] **Step 2: Compile and score**

```python
# Add to main():
body = author_small_talk()
path = OUT_DIR / "SmallTalk_AhEe_v1.json"
path.parent.mkdir(parents=True, exist_ok=True)
with open(path, "w") as f:
    f.write(body.to_compiled_json(provenance="shipping-forge"))
print(f"Saved: {path}")
```

Run: `python -m pyruntime.recipes.author_shipping`

- [ ] **Step 3: Validate with trench-mcp**

Run: `mcp__trench-mcp__analyze_state(cartridge_path="cartridges/candidates/SmallTalk_AhEe_v1.json", morph=0.5, q=0.5)`

**Acceptance criteria:**
- talkingness > 0.5
- ridge prominence > 0.5
- formant separation > 0.5 octaves
- vowel bias matches "ah" or "eh" at midpoint
- No reactor stress above 0.5

- [ ] **Step 4: Score with trajectory fitness**

Run: `python -c "from pyruntime.forge_optimize import _eval_body; ..."`

**Acceptance criteria:**
- trajectory > 0.3
- continuity > 0.5
- occupancy > 0.6
- gate_fail = False

- [ ] **Step 5: Iterate radii/zero offsets if criteria not met**

Adjust `radius` values (0.985-0.998 range) and `zero_offset_octaves` (0.5-2.0 range) until both MCP and trajectory scores pass. The reference body Ooh_to_Eee_(approx) has r=0.9917-0.9995.

- [ ] **Step 6: Commit passing body**

```bash
git add pyruntime/recipes/author_shipping.py cartridges/candidates/SmallTalk_AhEe_v1.json
git commit -m "feat: Small Talk Ah-Ee v1 — vocal morph body passing trajectory/MCP gates"
```

---

### Task 3: Author Speaker Knockerz (body 1)

**Acoustic targets (from SHIPPING.md):**
- Sub anchor below 60 Hz that NEVER disappears
- Morph path: Vault → Chest Resonance → Choke → Cardboard Rip → Rattle → Cone Cry → Total Fracture
- Invariant: fundamental below 60 Hz stays phase-locked across entire morph

**Construction pattern:**
- Stage 0: Sub anchor at 40-50 Hz (r=0.998, PURE contour — no zero, just pole)
- Stage 1: Low-mid resonance 150-400 Hz (r=0.990, INTERIOR_ZERO)
- Stage 2: Mid-band stress 400-1200 Hz (r=0.985, UNIT_CIRCLE zero for notch)
- Q axis: sub anchor holds (same r at Q0 and Q100), upper stages tighten

- [ ] **Step 1: Implement author_speaker_knockerz()**

Key design: M0 = "Vault" (sub + clamped highs), M100 = "Total Fracture" (sub + violent HF).
Q0 = relaxed, Q100 = pressured (sharper mid resonances but sub unchanged).
The sub anchor stage (S0) has identical params across all 4 corners.

- [ ] **Step 2: Compile and validate**

**Acceptance criteria:**
- Sub energy (20-60 Hz) stays within 3 dB across all 25 grid points (5x5)
- morph_distance > 5.0
- ridge > 0.5 at morph endpoints

- [ ] **Step 3: Validate with trench-mcp analyze_state at M0, M50, M100**

- [ ] **Step 4: Commit**

---

### Task 4: Author Aluminum Siding (body 2)

**Acoustic targets (from SHIPPING.md):**
- Permanent 1 kHz void (midrange scooped at all morph/Q positions)
- All energy focused in extreme highs (7-18 kHz)
- Morph path: Dull Silver → Sheen → Glass Stress → Sibilant Fold → Aluminium Tear → Dog Whistle → Shatter Point

**Construction pattern:**
- Stage 0: High shelf/peak at 10-12 kHz (r=0.975, UNIT_CIRCLE zero at ~1 kHz for the void)
- Stage 1: Moving resonance 7-18 kHz (r=0.980-0.990)
- Stage 2: Anti-resonance notch locked at 1 kHz (r=0.995, zero_forced at 1 kHz)
- Q axis: high-band peaks sharpen, void deepens

- [ ] **Step 1: Implement author_aluminum_siding()**

- [ ] **Step 2: Compile and validate**

**Acceptance criteria:**
- 1 kHz region below mean at all grid points (the void)
- Peak energy above 5 kHz at morph > 0.5
- morph_distance > 5.0

- [ ] **Step 3: Validate with trench-mcp**

- [ ] **Step 4: Commit**

---

### Task 5: Author Cul-De-Sac (body 4)

**Acoustic targets (from SHIPPING.md):**
- M0 = thick tube resonance (single dominant peak)
- M50 = State Boundary / phase cancellation null
- M100 = scattered comb filter (multiple narrow peaks)
- Invariant: root hum at fundamental stays present

**Construction pattern (most complex body):**
- M0 corners: 1-2 active stages, concentrated energy (thick resonance)
- M100 corners: 4-5 active stages, spread across spectrum (comb)
- The bilinear interpolation between few-stage and many-stage creates the "fracture" at midpoint
- Q axis: tube narrows (higher Q on dominant peak) vs comb spreads

- [ ] **Step 1: Implement author_cul_de_sac()**

Key insight: stages that are passthrough at M0 and active at M100 naturally create the "fracture" transition through bilinear interpolation. The coefficients blend from passthrough (c0=1,c1-c4=0) to resonator, creating emergent peaks at intermediate morph positions.

- [ ] **Step 2: Compile and validate**

**Acceptance criteria:**
- Dynamic range at M0 < 20 dB (concentrated)
- Dynamic range at M100 > 30 dB (scattered)
- Root frequency (20-200 Hz) stays present across all grid points

- [ ] **Step 3: Validate with trench-mcp**

- [ ] **Step 4: Commit**

---

### Task 6: Final scan and ranking

- [ ] **Step 1: Run scan_vault on all 4 candidates**

Run: `mcp__trench-mcp__scan_vault(vault_dir="cartridges/candidates")`

- [ ] **Step 2: Compare against P2K reference bodies**

All 4 shipping bodies should score composite > 0.60 (minimum threshold).
Small Talk Ah-Ee should target > 0.80 (matching Early Rizer/Ooh_to_Eee tier).

- [ ] **Step 3: Hotswap each body for live audition**

Run: `mcp__trench-mcp__hotswap(body_path="cartridges/candidates/SmallTalk_AhEe_v1.json")`

Repeat for each body. Listen through the plugin.

- [ ] **Step 4: Promote passing bodies to cartridges/ root**

```bash
cp cartridges/candidates/SmallTalk_AhEe_v1.json cartridges/SmallTalk_AhEe.json
# repeat for each passing body
```

- [ ] **Step 5: Final commit**

```bash
git add cartridges/*.json
git commit -m "feat: 4 shipping bodies — Speaker Knockerz, Aluminum Siding, Small Talk Ah-Ee, Cul-De-Sac"
```

---

### Known risks

1. **Radius tuning is empirical.** P2K reference bodies use r=0.975-0.999 but the exact values depend on the frequency band. Lower frequencies need higher radii for equivalent perceptual sharpness. Plan for 2-3 iterations per body.

2. **Zero placement affects gain structure.** Interior zeros can attenuate or boost the cascade depending on their position relative to the pole. The `zero_offset_octaves` parameter needs tuning per stage.

3. **The "fracture" transition in Cul-De-Sac** depends on how passthrough→resonator blending behaves. Test intermediate morph positions (0.4-0.6) carefully for the State Boundary null.

4. **Heritage integer templates are NOT usable as-is** for shipping bodies due to the gain collapse in `_fw_words_to_kernel`. They provide design intent (stage patterns, frequency ranges) but not shipping coefficients.
