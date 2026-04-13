# TRENCH SESSION STATE — 2026-04-13

Context reset dump. Full content, no compaction. Read this before resuming.

## Workspace map (three directories, three domains)

| Path | Role | Git |
|---|---|---|
| `C:/Users/hooki/Trench` | Domain 1 — active TRENCH implementation (this session's focus) | root git (branch `codex/v1`) |
| `C:/Users/hooki/trenchwork_clean` | Prior reference implementation with full `trench-core::engine::FilterEngine` + AGC + DC blocker | own git |
| `C:/Users/hooki/trench_re_vault` | Domain 2 — firmware RE data, ghidra dumps, captured reference wavs, stage-response models | own git |
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
- `cartridge.schema.json` — `compiled-v1` wire format, validates 46 live cartridges
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
- Validator: `cartridge.schema.json`. 46 cartridges in `cartridges/` validate.
- Two parallel layouts exist in the wild:
    - Wire (keyframe-object form): `{format, name, sampleRate, keyframes: [{label, morph, q, boost, stages: [{c0..c4}]}]}`. Trench cartridges use this.
    - Array form (legacy): `{version, corners: {M0_Q0: [[c0..c4], ...], ...}}`. `trench-core/src/cartridge.rs::Cartridge::from_json` also accepts this; schema does not.

**Hard bans (enforcement, not suggestion):**
- No RBJ cookbook for character filters.
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

**Architectural truth from DillusionMan 2005 Peak/Shelf Morph tutorial
(cross-checked this session):** The E-mu Peak/Shelf Morph filter does
NOT use a fixed 7th global lowpass stage. The lowpass is authored into
the same 6 morph stages via the `SHELF` parameter (-64 = LP, +63 = HP).
At M0 the top stages are written as an LP; at M100 they warp to
HP/shelf territory. This validates the existing TRENCH 6-active-stage
architecture — no 7th stage needed. Use this when authoring Speaker
Knockerz / Cul-De-Sac bass behavior: encode the LP in the M0 corner of
stages 1–2, explode to shelf/HP at M100. Maps to heritage P2k_021.

## Current work state (uncommitted)

**`Trench/` (root git, branch `codex/v1`):**
- Modified: CLAUDE.md, SPEC.md, README.md, AGENTS.md, handoff.json, docs/body_authoring_seed.md, docs/world_body_candidate_schema.md, docs/trench_master_prompt_library.md
- Added: DOCTRINE.md, MODES.md, PHONEMES.md, cartridge.schema.json, check, tools/parity_null.py, trench-core/tests/null_targets.rs, docs/archive/2026-04-05-trench-v1-instrument-architecture-design.md, docs/archive/qsound_spatial.md, this file
- Renamed: SHIPPING.md → BODIES.md
- Deleted: REPO_TRUTH.json, REPO_TRUTH.md, tools/update_repo_truth.py, docs/sonic_tables/{heritage_phonemes.json, sonic_tables_unified_bark.json, tables_bark.json}
- Many untracked new dirs (docs/calibration, docs/identity, docs/looperator, vault/_atlas, etc.) from prior sessions — left alone.

**`Trench/trench-juce/plugin/.git` (plugin subtree, branch `codex/monorepo-engine`):**
- Last landed commit `b71c854 sync plugin workspace: engine subtree, gui/thumbwheel updates, sync tools`
- Staged deletion: 197 files in `engine/` subtree (was a duplicate of root `trench-core/`, confirmed identical source trees, no build wiring referenced it)
- `.gitignore` updated to exclude `*.blend1`

**`Trench/trench-juce/.git`**: deleted this session. It was a dormant empty repo with zero commits. Plugin's own `.git` at `trench-juce/plugin/.git` is intact.

## Session focus

User goal was a "hyper optimised lean workspace" followed by real hardware
parity measurement against firmware RE reference captures. Everything
below tracks that focus.

### Verification: `./check` pipeline (all passing)

1. Canonical doc set present
2. No orphan `SHIPPING.md`
3. `cargo check --workspace --tests` (type check, no link — avoids Git-Bash `/usr/bin/link` shadowing MSVC `link.exe`; check script prepends MSVC to PATH)
4. `pyruntime` imports
5. Passthrough null targets (rust test `trench-core/tests/null_targets.rs`): static and dynamic both at -300 dBFS (bit-exact passthrough)
6. Hardware parity null (report-only, `tools/parity_null.py`)
7. 46 compiled-v1 cartridges validate against `cartridge.schema.json`

### Hardware parity null — current best numbers

Pipeline: `raw-ROM stage → SOS cascade (scipy.signal.sosfilt) → AGC (ported) → ×boost → lag+gain null vs hardware capture`

Inputs:
- Cartridge: `trenchwork_clean/cartridges/00_talking_hedz.json` (raw-ROM stage layout with `{a1, r, val1, val2, val3, flag}` per stage, `name: P2k_013`, `sampleRate: 39062.5` field but coefficients are baked for 44100 via `type3_to_encoded` default)
- Dry input: `trenchwork_clean/ref/bypassed-pinknoise.wav` (132694 samples, 44100 Hz, pink noise through bypassed filter)
- References: `trenchwork_clean/ref/hedz*.wav` (same session as the dry)

Results:

| corner     | lag | best-fit gain | null rms abs | ref rms   | rel null |
|------------|----:|--------------:|-------------:|----------:|---------:|
| M0_Q0      |   8 | 0.5013        | -108.9 dBFS  | -32.1 dBFS| **-76.8 dB** |
| M0_Q100    |   0 | 0.9989        |  -51.8 dBFS  | -19.8 dBFS| **-32.0 dB** |
| M100_Q0    |   0 | 1.0002        |  -59.5 dBFS  | -25.9 dBFS| **-33.6 dB** |
| M100_Q100  |   0 | 0.9962        |  -44.2 dBFS  | -19.3 dBFS|  -24.9 dB    |

M0_Q0 nulls shape-exact. Other three null shape-close but not clean.

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

### AGC algorithm (ported from `trenchwork_clean/trench-core/src/agc.rs`)

Documented in the source as verified against the EmulatorX.dll binary.
16-entry table, per-sample execution, gain never exceeds 1.0, no floor.

```
AGC_TABLE = [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160,
             0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120]

def agc_step(sample, agc_gain):
    idx = int(agc_gain * abs(sample)) & 0xF
    new_gain = agc_gain * AGC_TABLE[idx]
    agc_gain = new_gain if new_gain < 1.0 else 1.0
    return sample * agc_gain, agc_gain
```

Applying this to the raw cascade output closes the null by 23–63 dB
depending on corner. The 2026-03-12 RE note
(`talking_hedz_x3_surface_probe_20260312.py`) had already flagged the
gap as "gain staging missing" — AGC is that missing piece.

### What is still wrong (analysis)

The remaining gap splits into two questions:

**A. Why does M0_Q0 null shape-perfectly at a 2.0× amplitude mismatch
(best-fit gain 0.5013)?**

M0_Q0 is the quietest corner (ref rms -32 dBFS, lowest Q, lowest morph).
AGC barely touches quiet signals (index 0/1 of the table are 1.0001),
so AGC isn't the source of the 2× factor. Possible causes, in order of
likelihood:

1. **Untracked per-corner pre-AGC scale.** `FilterEngine` does more than
   cascade+AGC+boost; it has motion-pursuit smoothing between desired
   and actual coefficients (`alpha` per stage) and a coefficient ramp
   from passthrough to target on initial load. For a static corner with
   no prior state, the first 32 samples differ from steady-state sosfilt
   output. However this is a transient, not a global 2× — unlikely the
   cause.
2. **A control-rate envelope in the cartridge that the raw ROM json
   carries separately.** The wire format has `boost` per keyframe which
   I apply; there may be another per-corner scalar (e.g. the first-block
   pursuit pushes the desired coefficients toward target over the block
   while the state is still in the passthrough starting condition,
   producing a time-varying gain that a scalar best-fit averages out to
   0.5).
3. **Reference hedzm0q0.wav is shorter (101280 samples vs 132694 for
   the other three).** Different trim could imply different recording
   session or different processing — a 0.5× gain normalization in the
   capture chain would match exactly. Worth verifying against capture
   manifest.

**B. Why do M0_Q100, M100_Q0, M100_Q100 null at -25 to -33 dB instead of
closing?**

Gain is ~1.0 on those three, so amplitude is correct but shape has a
residual error. Candidates:

1. **Coefficient ramping.** Real engine ramps coefficients from the
   starting state (passthrough) to the corner target over 32 samples.
   sosfilt is steady-state from sample 0. First block of prediction vs
   reference will differ. For 132694 samples this is small fractional
   energy but could set a floor near -40 dB.
2. **Motion pursuit alpha filter.** engine.rs processes
   `desired → motion` via alpha smoothing per stage. My sosfilt uses
   the desired coefficients directly. On a static corner this
   converges quickly, but the convergence trajectory is not identical.
3. **Minifloat quantization.** The wire path goes through packed u16
   minifloats (see `EncodedCorners.packed_m*`). My cascade uses the
   unquantized `val1/val2/val3` floats. The hardware output was
   rendered after minifloat round-trip, producing quantization noise
   present in the reference but not my prediction.
4. **DC blocker in reference.** The hardware reference may have been
   rendered with DC blocker enabled; my sosfilt output has none. I
   tried adding my own DC blocker — broke alignment (phase shift).
   Correct fix: either both sides get identical DC blocker, or both
   skip it.

**Previous session notes:** The user states "we null tested previously"
— implying a prior session produced a clean null. The
`trench_re_vault/scratch/talking_hedz_null_test_2026-03-14/` directory
has an `src/main.rs` that uses `trench_core::engine::FilterEngine` with
two modes (`Production` with AGC + clamp, `Linear` with both disabled).
It writes reports to
`scratch/talking_hedz_null_test_2026-03-14/reports/talking_hedz_null_test_report.md`.
That binary produces the real null numbers. It has not been run from
this session; its report file was not inspected. **Running it is the
single next step most likely to reveal what my python pipeline is
missing.**

## Blockers

1. `Trench/trench-core/` has no `engine.rs`, no `agc.rs`, no dc-blocker —
   the FilterEngine was never ported from `trenchwork_clean`. Current
   parity null runs entirely in python; the rust workspace cannot
   reproduce it.
2. M0_Q0 2× amplitude mystery is unresolved. Until it is, even a perfect
   port won't explain the 0.5 best-fit gain.
3. Staged work is uncommitted on both Trench root and plugin repo. A
   context reset risks losing the mental map of what was staged vs
   working. This file is the recovery map.
4. `tools/parity_null.py` depends on `trenchwork_clean/ref` and
   `trenchwork_clean/cartridges/00_talking_hedz.json`. It soft-skips if
   those paths are missing — CI on a clean machine will not run it.
5. Hardware parity step in `./check` is report-only. It prints numbers
   but does not fail on regression. No baseline file to compare against.

## Next steps (ordered)

1. **Run the existing null_test binary.** Build and run
   `trench_re_vault/scratch/talking_hedz_null_test_2026-03-14/experiments/talking_hedz_null_test`.
   Read `reports/talking_hedz_null_test_report.md`. The numbers there
   are the ground truth for what a working null should look like.
   Compare corner by corner against the python numbers in this file.
   The delta isolates exactly what the python pipeline is missing
   (motion pursuit? minifloat quantization? DC blocker state?).
2. **Inspect `hedzm0q0.wav` manifest.** Find `CAPTURE_MANIFEST.json` in
   `trench_re_vault/datasets/talking_hedz_bracketed_residual/2026-03-12`
   and determine whether M0_Q0 was recorded with different gain,
   different trim, or different session. If yes, the 0.5 best-fit gain
   is a capture-side artifact and can be documented rather than chased.
3. **Port FilterEngine + AGC + DC blocker** from `trenchwork_clean/trench-core/src/{engine.rs, agc.rs}` into `Trench/trench-core/src/`. Keep the Trench cascade
   frozen per doctrine; the port is a new `engine` module on top of
   existing `Cascade`. Add a rust integration test that reproduces the
   python parity numbers, wire it into `./check` replacing the
   report-only python step.
4. **Commit the staged work** on both repos. Two commits:
   - Trench root: lean workspace + null pipeline + session state
   - plugin: delete `engine/` duplicate subtree
5. **Audition the DillusionMan LP-in-morph finding.** Author a small
   candidate body (stage 1–2 as LP at M0, shelf at M100) and render
   through `./check`'s cartridge schema + (eventually) the ported
   engine. Compare spectrally to heritage P2k_021. Confirms the "no
   7th stage" architectural commitment holds for Speaker Knockerz bass
   behavior without a separate global LP.

## Operator handoff cue

When the operator resumes, start with step 1 above. Do not re-derive
the AGC; do not re-search for reference files; do not re-explore the
workspace map. All of that is above.
