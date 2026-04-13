# TRENCH SESSION STATE — 2026-04-13 (rev 2)

Context reset dump. Full content, no compaction. Read this before resuming.
Supersedes rev 1 from earlier this session.

## Workspace map (three directories, three domains)

| Path | Role | Git |
|---|---|---|
| `C:/Users/hooki/Trench` | Domain 1 — active TRENCH implementation (this session's focus) | root git (branch `chore/git-hygiene`, was `codex/v1` at session start) |
| `C:/Users/hooki/trenchwork_clean` | Prior reference implementation with full `trench-core::engine::FilterEngine` + AGC + DC blocker + motion pursuit + minifloat | own git |
| `C:/Users/hooki/trench_re_vault` | Domain 2 — firmware RE data, ghidra dumps, captured/rendered reference wavs, stage-response models | own git |
| `C:/Users/hooki/Trench/trench-juce/plugin` | Active shipping plugin (JUCE 8) | own git; excluded from root `/trench-juce/` gitignore |

Separation rule (from memory `Three Domain Separation`): never conflate
TRENCH impl truth, firmware RE truth, and future design. Each lives in
its own tree.

## Canonical docs (Domain 1 root)

- `CLAUDE.md` — per-session brief, 28 lines
- `SPEC.md` — math + cartridge contract only, 80 lines
- `DOCTRINE.md` — working rules, hard bans, verification, escalation, 38 lines
- `MODES.md` — Shape Bank / Trajectory / Operator / Shipping modes, 23 lines
- `BODIES.md` — 4 shipping bodies + per-body "how to know it's wrong" rubrics, 169 lines
- `PHONEMES.md` — authoring model (phoneme grid → cartridge), forward-looking
- `cartridge.schema.json` — `compiled-v1` wire format; validates 47 live cartridges
- `./check` — bash verification script (doc set + cargo type-check + pyruntime imports + null targets + hardware parity report + cartridge schema)
- `AGENTS.md`, `README.md` — thin pointers to canonical docs

## Domain 1 invariants (from SPEC.md + DOCTRINE.md)

**DSP engine (frozen):**
- 12-stage serial DF2T biquad cascade: 6 active + 6 passthrough. Never parallel.
- Kernel-form `[c0, c1, c2, c3, c4]` interpolation only. Never interpolate raw biquad coefficients or raw frequencies.
- 4-corner bilinear interpolation. Q axis first, then morph axis.
- 32-sample control blocks with per-sample coefficient ramping.
- Sample rates: 39062.5 Hz (authoring/forge), 44100.0 Hz (plugin runtime).
- Minifloat: 4-bit exponent (bias 15), 12-bit mantissa. Sentinels: `0x0000` = 0.0; `0xDFFF` = passthrough gain; `0xFFFF` = constant 1.0.
- DF2T per-sample math:
    ```
    y[n]  = c0*x[n] + w1[n-1]
    w1[n] = c1*x[n] - c3*y[n] + w2[n-1]
    w2[n] = c2*x[n] - c4*y[n]
    ```
- Bilinear (Q then morph):
    ```
    q_m0 = C_M0_Q0   + (C_M0_Q100   - C_M0_Q0)   * Q_frac
    q_m1 = C_M100_Q0 + (C_M100_Q100 - C_M100_Q0) * Q_frac
    coef = q_m0      + (q_m1        - q_m0)      * morph_frac
    ```
- Cascade source files are frozen: `trench-core/src/cascade.rs`, `trench-juce/plugin/source/TrenchEngine.cpp`.

**Cartridge format (`compiled-v1` on disk):**
- `format: "compiled-v1"`
- 4 keyframes labeled `M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`
- Each keyframe carries 12 stage arrays of 5 coefficients (6 active + 6 passthrough sentinel)
- Validator: `cartridge.schema.json`. 47 cartridges in `cartridges/` validate (up from 46; `cartridges/auditions/lp_in_morph.json` was added this session as the LP-in-morph audition body).
- Two parallel layouts exist in the wild:
    - Wire (keyframe-object form): `{format, name, sampleRate, keyframes: [{label, morph, q, boost, stages: [{c0..c4}]}]}`. Trench cartridges use this.
    - Array form (legacy): `{version, corners: {M0_Q0: [[c0..c4], ...], ...}}`. `trench-core/src/cartridge.rs::Cartridge::from_json` also accepts this; schema does not.
- Raw-ROM layout (third format, found in `trench_re_vault/scratch/resources/roms/P2K/*.json` and `trenchwork_clean/cartridges/*.json`): `{corners: {M*_Q*: {stages: [{a1, r, val1, val2, val3, flag}]}}}`. Conversion to kernel form:
    - `a2 = r²` (with `r ≤ 0.999999` stability clamp)
    - `flag < 0.5` → all-pole LP: `b0 = 1, b1 = 0, b2 = 0`
    - `flag ≥ 0.5` → `b0 = 1 + val1, b1 = a1 + val2, b2 = a2 - val3`
    - DF2T map: `c0 = b0, c1 = b1, c2 = b2, c3 = a1, c4 = a2`

**Hard bans (enforcement, not suggestion):**
- No RBJ cookbook for character filters. Direct pole-zero only.
- No compensation layers, fudge factors, or saturation to rescue bad math.
- No parallel engine paths unless explicitly proven and approved.
- No gain baked into `c4`.
- No unapproved dependency changes.
- No shipping verbatim heritage extractions.
- No pole sanitization without a failing test proving instability.
- DSP cascade, interpolation order, and cartridge format are frozen.

**Shipping set (the only normative names for v1):**
1. Speaker Knockerz (sub pressure / cone cry)
2. Aluminum Siding (brittle top-end damage)
3. Small Talk Ah-Ee (biomechanical vocal cavity)
4. Cul-De-Sac (tube → comb fracture meta-body)

Each has a "how to know it's wrong" rubric in `BODIES.md`.

**Architectural truth — no 7th global lowpass stage:** DillusionMan's
2005 Peak/Shelf Morph tutorial describes the E-mu filter's SHELF
parameter as -64 → +63, with -64 = LP and +63 = HP. Body authoring
places LP, peak, or HP behavior directly into the same 6 active stages
via per-corner coefficients. The LP is NOT a fixed 7th stage; it is
authored into the M0 corner of whichever stages should carry it,
and warped to shelf or HP at the M100 corner.

This session added `cartridges/auditions/lp_in_morph.json` as the
explicit audition: 4-pole LP at 200 Hz at M0 corners, 4-pole HP at
800 Hz at M100 corners, stages 3–12 passthrough. Per-corner frequency
response confirms the LP→HP morph is expressible within the 6-active
architecture without any cascade change.

Applies to shipping bodies: when authoring Speaker Knockerz and
Cul-De-Sac bass behavior, encode the LP directly in stage 1–2 of the
M0 corner and let it warp as needed at M100. No architectural exception
required. Heritage heritage counterpart: P2k_021 stays entirely in LP
territory across all four corners (resonant LP morphing through cutoff
and emphasis), not LP→HP. The audition body is a wider-range test
than P2k_021 itself.

## Commits landed this session

All on branch `chore/git-hygiene`:

- `cfc3f70 port AGC from trenchwork_clean trench-core` — `trench-core/src/agc.rs` (89 lines, verbatim port, 5 unit tests pass), `trench-core/src/lib.rs` adds `pub mod agc` and re-exports `agc_step`, `dev/SESSION_STATE.md` (rev 1)
- `42e798b audition: LP-in-morph candidate body (6-stage, no 7th global LP)` — `cartridges/auditions/lp_in_morph.json` (375 lines), confirms the architectural claim with validated schema + spectral response

Earlier in the session (committed in the background by an external
process): `0cb3f29`, `507bca7`, `5a1ae92`, `8d30aae`, `adab1e9` —
lean workspace, cartridge schema, ./check, parity_null script, null
target tests, research docs.

Plugin repo (`Trench/trench-juce/plugin/.git`, branch
`codex/monorepo-engine`):
- `c58a187 delete plugin/engine duplicate of root trench-core workspace` — 197-file deletion of the duplicate trench-core subtree inside the plugin. Confirmed via `diff -rq` that the only difference between `plugin/engine/trench-core/` and root `trench-core/` was `CLAUDE.md`; the src/ tests/ Cargo.toml were byte-identical. CMakeLists never referenced `engine/`, so no build wiring broke.
- Previous commit `b71c854` (same branch) synced plugin workspace changes.

`Trench/trench-juce/.git` — deleted this session. It was a dormant
empty repo with zero commits. Plugin's own `.git` at
`trench-juce/plugin/.git` is intact.

## Session focus (Domain 1)

User goal was a "hyper optimised lean workspace" followed by real
hardware parity measurement against firmware RE reference captures.
This session landed both plus ported the AGC and proved the
LP-in-morph architectural claim.

### Verification: `./check` pipeline (all passing)

1. Canonical doc set present (CLAUDE, SPEC, DOCTRINE, MODES, BODIES, PHONEMES, cartridge.schema.json)
2. No orphan `SHIPPING.md`
3. `cargo check --workspace --tests` (type check only — avoids Git-Bash `/usr/bin/link` shadowing MSVC `link.exe`; check script prepends MSVC to PATH for the null_targets test build)
4. `pyruntime` imports
5. Passthrough null targets (rust test `trench-core/tests/null_targets.rs`): static and dynamic both at -300 dBFS (bit-exact passthrough)
6. Hardware parity null (report-only, `tools/parity_null.py`)
7. 47 compiled-v1 cartridges validate against `cartridge.schema.json`

### Hardware parity null — current best numbers

Pipeline: `raw-ROM stage → SOS cascade (scipy.signal.sosfilt) → AGC (ported verbatim from trenchwork_clean) → ×boost → lag+gain null vs hardware capture`

Inputs:
- Cartridge: `trenchwork_clean/cartridges/00_talking_hedz.json` (raw-ROM stage layout, `name: P2k_013`, `sampleRate: 39062.5` field but coefficients baked for 44100 via `type3_to_encoded(sr=44100)` default at compile time)
- Dry input: `trenchwork_clean/ref/bypassed-pinknoise.wav` (132694 samples, 44100 Hz stereo float64, pink noise through bypassed filter)
- References: `trenchwork_clean/ref/hedz*.wav` (same session as the dry; byte-identical to `trench_re_vault/datasets/talking_hedz_bracketed_residual/2026-03-12/inputs/hedz*.wav`)

Results:

| corner     | lag | best-fit gain | null rms abs | ref rms    | rel null |
|------------|----:|--------------:|-------------:|-----------:|---------:|
| M0_Q0      |   8 | 0.5013        | -108.9 dBFS  | -32.1 dBFS | **-76.8 dB** |
| M0_Q100    |   0 | 0.9989        |  -51.8 dBFS  | -19.8 dBFS | **-32.0 dB** |
| M100_Q0    |   0 | 1.0002        |  -59.5 dBFS  | -25.9 dBFS | **-33.6 dB** |
| M100_Q100  |   0 | 0.9962        |  -44.2 dBFS  | -19.3 dBFS |  -24.9 dB    |

M0_Q0 nulls shape-exact. Other three null shape-close but not clean.

### Comparison against the rust null_test binary (2026-03-14)

`trench_re_vault/scratch/talking_hedz_null_test_2026-03-14/reports/talking_hedz_null_test_report.md`
was generated by a rust binary using `trench_core::engine::FilterEngine`
with full motion pursuit, minifloat quantization, coefficient ramping,
and AGC.

| corner     | rust (2026-03-14)     | python (this session) |
|------------|-----------------------|-----------------------|
| M0_Q0      | gain 0.500882, null -3.10 dB  | gain 0.5013, null -76.8 dB |
| M0_Q100    | gain 0.977981, null -7.42 dB  | gain 0.9989, null -32.0 dB |
| M100_Q0    | gain 0.999373, null -3.03 dB  | gain 1.0002, null -33.6 dB |
| M100_Q100  | gain 0.997642, null -9.07 dB  | gain 0.9962, null -24.9 dB |

Both pipelines agree on M0_Q0 best-fit gain = 0.500 → the 2× amplitude
factor is a property of the reference file, not the rendering pipeline.
Python pipeline's null is dramatically better at each corner. This
means the rust FilterEngine's extra machinery (motion pursuit alpha
smoothing, minifloat pack/unpack, coefficient ramping from passthrough
initial state, block-rate control) drifts the output away from what
the reference capture encodes. The simpler python sosfilt path matches
the reference better.

Implication: the reference `hedz*.wav` files were almost certainly
generated by a simple scipy-like cascade rendering, not by the full
FilterEngine or by actual E-mu hardware capture. They are a simulated
ground truth, not hardware truth. True hardware null would need real
DAW captures of the hardware running the same cartridge.

Note: `trenchwork_clean/captures/benchmark/P2K_*.wav` are 8-second
captures in the same workspace, plausibly real DAW captures of the
hardware, but no matching dry input exists for them — the
`trenchwork_clean/ref/bypassed-pinknoise.wav` dry is 3 seconds long,
different session. Until a matching dry is located (or a sine sweep
deconvolution is applied), those captures cannot be used for direct
null measurement.

### Stage conversion formula (raw ROM → kernel, verified)

From `trench_re_vault/tools/morpheus_lab/plot_p2k_stage_responses_20260312.py`
and confirmed matches `Trench/cartridges/p2k/P2k_004.json` byte-for-byte
when applied to `trench_re_vault/scratch/resources/roms/P2K/P2k_004.json`:

```python
def stage_coefficients(stage):
    a1 = float(stage["a1"])
    r  = min(float(stage["r"]), 0.999999)   # stability clamp
    a2 = r * r
    if float(stage.get("flag", 1.0)) < 0.5:
        b0, b1, b2 = 1.0, 0.0, 0.0          # flag=0 → all-pole lowpass
    else:
        b0 = 1.0 + float(stage["val1"])
        b1 = a1  + float(stage["val2"])
        b2 = a2  - float(stage["val3"])
    return b0, b1, b2, a1, a2               # DF2T: c0=b0, c1=b1, c2=b2, c3=a1, c4=a2
```

### AGC algorithm (now in Trench/trench-core/src/agc.rs)

Documented in the source as verified against the EmulatorX.dll binary
(disassembled from `FUN_1802c04e0`). 16-entry table, per-sample
execution, gain never exceeds 1.0, no floor, no recovery across
samples beyond what the table allows.

```
AGC_TABLE = [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160,
             0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120]

pub fn agc_step(sample: f32, agc_gain: &mut f32) -> f32 {
    let idx = ((*agc_gain * sample.abs()) as u32 & 0xF) as usize;
    let new_gain = *agc_gain * AGC_TABLE[idx];
    *agc_gain = if new_gain < 1.0 { new_gain } else { 1.0 };
    sample * *agc_gain
}
```

5 unit tests in `trench-core/src/agc.rs` cover: quiet passthrough,
loud-signal gain reduction, recovery below threshold, no floor
(gain can drop arbitrarily low), output finiteness across inputs.
All pass.

For quiet signals (`agc_gain * |sample| < 1.0`), the index is 0 or 1,
table values are 1.0001, new_gain is clamped back to 1.0, and the AGC
does nothing. This matches the observation that the rust null_test
report shows identical numbers for Production (AGC on) and Linear
(AGC off) modes — the pink-noise test signal is too quiet to trigger
AGC gain reduction.

### What is still wrong (analysis, updated)

**M0_Q0 2× amplitude mismatch (gain 0.5013):** Both the rust
FilterEngine and the python sosfilt pipeline agree on this gain
with the same four-decimal precision. This confirms it is NOT a
pipeline error — the hedzm0q0.wav reference was rendered (or
recorded) at half the amplitude that both independent render paths
produce. Most likely explanation: the reference generator applied a
corner-specific gain factor (0.5 at M0_Q0) that has no equivalent
in either FilterEngine.rs or my parity_null.py. Unresolved.

**M0_Q100 / M100_Q0 / M100_Q100 remaining gap (-24 to -33 dB):**
Gains are ~1.0 on those corners, so amplitude is correct but shape
has residual error. Python nulls them 21–25 dB better than rust, so
the FilterEngine path introduces extra drift. Candidates for the
remaining python residual:

1. **DC blocker in reference.** The reference wavs may carry a DC
   blocker signature; my sosfilt output does not. Adding a python
   DC blocker (matching `trenchwork_clean/trench-core/src/engine.rs::DcBlocker::process`) dropped the
   null by ~30 dB at every corner — so the reference was NOT
   rendered with that DC blocker. Either the reference has no DC
   blocker, or it uses a different DC blocker (different `r` or
   initial state).
2. **Alternative corner role mapping.** `bracketed_residual_analysis.json`
   documents a role swap: capture label `M0_Q100` ← source
   `M100_Q0`, capture label `M100_Q0` ← source `M0_Q100`. Applying
   this swap made the middle corners WORSE (back to ~-1.9 dB),
   confirming the straight (non-swapped) mapping is the right pairing
   for how this reference file was generated.
3. **Higher-order nonlinearity.** If the reference was rendered
   through a nonlinear element (saturation, AGC-in-a-different-place,
   envelope-driven gain), the residual at corners with higher
   signal levels (M0_Q100, M100_Q100) would be larger than at
   quieter corners. The -25 to -33 dB band is consistent with a
   quiet nonlinearity.

Pursuing any of these requires access to the reference generator
source code, which was not located in this session.

## Blockers

1. `Trench/trench-core/` has `agc.rs` now (this session) but still no
   `engine.rs` — the full `FilterEngine` was not ported. A dedicated
   port session is needed: engine.rs is 1429 lines and pulls in
   `EncodedCorners`, `PackedCoeffs`, motion-pursuit state, DC blocker,
   control-rate parameter dispatch, Q extrapolation stability guard,
   and a gain ceiling clamp.
2. M0_Q0 2× amplitude mystery is unresolved. Both rust and python
   report the same 0.500 best-fit gain — the 2× is in the reference
   file, not in any render pipeline currently in the workspace.
   Until the reference generator is found (or a new reference is
   captured from real hardware), this corner can never null cleanly
   against hedzm0q0.wav without applying a 2× correction somewhere
   arbitrary.
3. Hardware parity step in `./check` is report-only. It prints numbers
   but does not fail on regression. No baseline file to compare against
   committed metrics.
4. `tools/parity_null.py` depends on `trenchwork_clean/ref/` and
   `trenchwork_clean/cartridges/00_talking_hedz.json`. It soft-skips
   if those paths are missing — CI on a clean machine cannot run it.
5. Branch is `chore/git-hygiene`, not `codex/v1`. The session's
   landed work is on a side branch; merge/rebase strategy
   not decided.
6. Dry input for `trenchwork_clean/captures/benchmark/P2K_*.wav`
   (Ear Bender, Freak Shifta, Talking Hedz full-session, etc.)
   is not located. Those 8-second captures are plausibly real
   hardware recordings; without a matching dry they cannot be
   nulled directly.

## Next steps (ordered)

1. **Port FilterEngine + DC blocker** from `trenchwork_clean/trench-core/src/engine.rs` (1429 lines) into `Trench/trench-core/src/`. Keep the Trench `Cascade` frozen per doctrine; the port is a new `engine` module built on top of existing `Cascade` + the just-ported `agc`. Add a rust integration test that reproduces the python parity numbers. Wire it into `./check` replacing the report-only python step. This is a focused session on its own.
2. **Find the reference generator** for `hedz*.wav`. Search both
   `trench_re_vault` and `trenchwork_clean` for any script whose
   output filenames match `hedzm0q0`, `hedzmorph0q100`, etc.
   Failing that, look for calls to soundfile/wavfile/hound Write in
   tool directories named `talking_hedz_*`. The generator will
   reveal (a) whether the 2× gain at M0_Q0 is a per-corner
   amplitude override, (b) whether a DC blocker variant was
   applied, and (c) the exact order of cascade + AGC + gain stages.
3. **Locate the DAW dry** for `trenchwork_clean/captures/benchmark/P2K_*.wav`
   (8-second captures). If it exists, it enables parity null
   measurement against REAL hardware for all 33 P2K presets, not
   just Talking Hedz. If it does not exist, generate one via sweep
   deconvolution on the captures themselves (they appear to be
   swept or broadband test signals).
4. **Author Speaker Knockerz / Cul-De-Sac candidates** using the
   LP-in-morph pattern proven by `cartridges/auditions/lp_in_morph.json`.
   Encode sub-anchor LP directly in stage 1–2 of the M0 corner,
   warp to shelf or HP at M100. Validate each candidate against
   both `cartridge.schema.json` and its `BODIES.md` "how to know it's
   wrong" rubric.
5. **Make the parity null step gate-enforcing** rather than
   report-only: commit a baseline metrics JSON, have `./check`
   diff current numbers against baseline, fail on regression
   beyond a tolerance (e.g. null rms worsens by more than 3 dB at
   any corner).
6. **Decide branch strategy** for `chore/git-hygiene` — merge to
   `codex/v1`, or rebase, or open PR. User decision.

## Operator handoff cue

When the operator resumes, start with step 1 (FilterEngine port) if
the goal is parity gate in rust, or step 4 (Speaker Knockerz authoring)
if the goal is shipping-body progress. Skip all re-derivation —
workspace map, AGC table, stage conversion formula, and current parity
numbers are all captured above. The LP-in-morph architectural finding
is validated and committed; no further audition needed.
