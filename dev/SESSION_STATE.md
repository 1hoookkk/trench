# TRENCH SESSION STATE — 2026-04-13 (rev 3)

Context reset dump. Full content, no compaction. Read this before resuming.
Supersedes rev 2 from earlier this session.

Rev 3 adds: P2K source-of-truth resolution, Cheat Engine capture bug, canonical
reference regeneration, bit-exact parity gate covering all 35 skins.

## Workspace map (three directories, three domains)

| Path | Role | Git |
|---|---|---|
| `C:/Users/hooki/Trench` | Domain 1 — active TRENCH implementation (this session's focus) | root git (branch `chore/git-hygiene`) |
| `C:/Users/hooki/trenchwork_clean` | Prior reference implementation with full `trench-core::engine::FilterEngine` + AGC + DC blocker + motion pursuit + minifloat | own git |
| `C:/Users/hooki/trench_re_vault` | Domain 2 — firmware RE data, ghidra dumps, captured/rendered reference wavs, ROM binaries, stage-response models | own git |
| `C:/Users/hooki/Trench/trench-juce/plugin` | Active shipping plugin (JUCE 8) | own git; excluded from root `/trench-juce/` gitignore |
| `C:/Users/hooki/do-it/archive/data-dumps` | **Newly surfaced this session** — Cheat Engine dumps, `TalkingHedz*.json`, `zplane_kernel_raw.bin`, `ce_dump.txt`, `X3_TalkingHedz_Corners/` PNG screenshots. **Contains BROKEN CE stage captures that produced the legacy hedz*.wav references.** Not a ground-truth source. | separate tree |

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

**Trench compiled-v1 kernel mapping (confirmed this session via byte-level
agreement with raw ROM truth at 1e-7 precision for P2k_013):**
- `c0 = b0 = 1 + val1`
- `c1 = b1 = a1 + val2`
- `c2 = b2 = r² - val3`
- `c3 = a1` (standard polynomial coefficient, normally negative for stable LP)
- `c4 = a2 = r²`
- This is the **direct biquad form** of the DF2T math above, not the shifted
  minifloat-domain kernel form used by the older `trenchwork_clean`
  `pyruntime/encode.py::raw_to_encoded` (which stores `c0=2+b1/b0`, `c1=1-b2/b0`,
  `c2=2+a1`, `c3=1-r²`, `c4=b0`). **Two different conventions coexist in the
  workspace — do not confuse them.** Trench runtime uses the direct form.

**Cartridge format (`compiled-v1` on disk):**
- `format: "compiled-v1"`
- 4 keyframes labeled `M0_Q0`, `M0_Q100`, `M100_Q0`, `M100_Q100`
- Each keyframe carries 12 stage arrays of 5 coefficients (6 active + 6 passthrough sentinel)
- Validator: `cartridge.schema.json`. 47 cartridges in `cartridges/` validate.
- `Trench/cartridges/p2k/*.json` adds a top-level `provenance` field naming
  the source raw-ROM skin (e.g. `"P2k_013"`).
- Three parallel layouts exist in the wild:
    - **Wire compiled-v1 (keyframe-object form)**: `{format, name, provenance, sampleRate, stages, keyframes: [{label, morph, q, boost, stages: [{c0..c4}]}]}`. Trench p2k cartridges use this. Kernel mapping = direct biquad (above).
    - **Raw-ROM form**: `{name, filterType, stageCount, boost, corners: {M*_Q*: {stages: [{a1, r, val1, val2, val3, flag}]}}}`. Source of truth. Lives in `trenchwork_clean/datasets/p2k_skins/`, `trench_re_vault/scratch/resources/roms/P2K/`, and 3 other identical copies. Conversion to biquad:
        - `a2 = r²` (with `r ≤ 0.999999` stability clamp)
        - `flag < 0.5` → all-pole LP: `b0 = 1, b1 = 0, b2 = 0`
        - `flag ≥ 0.5` → `b0 = 1 + val1, b1 = a1 + val2, b2 = a2 - val3`
        - DF2T map: `c0 = b0, c1 = b1, c2 = b2, c3 = a1, c4 = a2`
    - **Array form (legacy)**: `{version, corners: {M0_Q0: [[c0..c4], ...], ...}}`. `trench-core/src/cartridge.rs::Cartridge::from_json` also accepts this; schema does not.

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

Audition body `cartridges/auditions/lp_in_morph.json` (4-pole LP at 200 Hz
at M0, 4-pole HP at 800 Hz at M100, stages 3–12 passthrough) confirms
this is expressible in 6 active stages without any cascade change.

## P2K source-of-truth chain (resolved this session)

The session-state rev 2 "M0_Q0 2× amplitude mystery" and "Q=100 bass residual"
were both artifacts of corrupted upstream stage data, not cascade math or
minifloat quantization.

### Canonical raw extraction (ground truth)

**sha256 prefix `9fb1bef0d0212980`** identifies the byte-identical canonical
raw extraction for `P2k_013.json` (4772 bytes, original extraction mtime
2026-02-17 08:52). Five byte-identical copies exist at:

1. `trench_re_vault/scratch/resources/roms/P2K/P2k_013.json`
2. `trench_re_vault/datasets/talking_hedz/2026-03-12/P2k_013.json`
3. `trench_re_vault/datasets/talking_hedz_bracketed_residual/2026-03-12/inputs/P2k_013.json`
4. `trenchwork_clean/datasets/p2k_skins/P2k_013.json`
5. (and copies in `trenchwork`, `trenchwork_backup`, `trenchwork_clean2`, `trenchwork_fresh`, `trenchwork_recovered`, plus 10+ `.claude/worktrees/agent-*/` clones)

Format: `{name, filterType, stageCount, boost, corners: {M0_Q0, M0_Q100, M100_Q0, M100_Q100: {stages: [{a1, r, val1, val2, val3, flag}]}}}`.

- **6 stages per corner** (not 5)
- `boost: 4.0` top-level (single corner-independent scalar)
- `flag=1` throughout for P2k_013 (all resonator, no all-pole)

`trenchwork_clean/datasets/p2k_skins/` is the canonical directory for 35
raw-ROM skins: 33 numbered (`P2k_000` through `P2k_032`) plus `Vocal_Ah_Ay_Ee`
and `Vocal_Oo_Ah`. Used as the source in `tools/render_canonical_refs.py`.

### The broken cartridge (confirmed corrupt)

`trenchwork_clean/cartridges/00_talking_hedz.json` (sha256 prefix
`cbb779d0f3558a08`, 5058 bytes, mtime 2026-03-23). Format has `keyframes`
list instead of `corners` dict.

**Corner label swap.** cart `M0_Q100` stage data == truth `M100_Q0`, and
cart `M100_Q0` stage data == truth `M0_Q100`. `M0_Q0` and `M100_Q100`
labels agree with truth. Verified by byte comparison of all 6 stages
of each corner.

Produced by an older `trenchwork_clean` compile path. Do not trust this
file. Do not load it in any parity pipeline.

### The Cheat Engine capture (also corrupt — separate bug)

`C:/Users/hooki/do-it/archive/data-dumps/TalkingHedz_Complete.json` is a
Cheat Engine / EmulatorX3 live memory capture from 2026-01-26 with
`captureDate` field. Three independent bugs vs canonical:

1. **Off-by-one stage indexing.** 5 stages instead of 6. Missed canonical
   stage 0 (the low-Q a1≈-0.14, r≈0.975 stage). `CE[i].a1 == canonical[i+1].a1`
   bit-exact for i ∈ 0..4.
2. **val1 captured as b0.** CE stored `1 + val1` (≈ 0.562) in the `val1`
   slot instead of raw val1 (≈ -0.438). `CE[i].val1 == 1 + canonical[i+1].val1`
   bit-exact.
3. **Radius drift.** `CE[i].r` differs from `canonical[i+1].r` by 0.003 to
   0.049. CE captured after some EmulatorX3 runtime decoding/interpolation
   step that perturbs radius. Canonical reads pre-decode from ROM.

Per-corner boost in CE ≈ 1.76 varying; canonical boost is flat 4.0.

**The legacy `trenchwork_clean/ref/hedz*.wav` references were rendered from
this broken CE capture**, not from canonical. That is the real reason the
rev 2 parity analysis hit a -42 dB floor at Q=100 corners. The cascade
math was correct; the references were wrong.

### The secondary truth chain (Morpheus, not P2K)

`trenchwork_clean/datasets/morpheus_zplane_library.json` — 289 z-plane cubes
decoded from `cubes_v1.01vc_170120.wav` (Rossum Morpheus z-plane library).
Each cube has 4 corners × 6 stages with `(r, freq_hz, a1, a2, r_byte, f_byte)`.
**Pure pole filters — no `val1/val2/val3` zero terms.** This is the Morpheus
base library, NOT the P2K extended set. P2K = Morpheus pole structure +
added zero placements via `val1/val2/val3`. Not a substitute for the P2K
raw extraction.

### The name table

`trenchwork_clean/datasets/p2k_filter_names.json` — 56 confirmed P2K filter
names extracted from EmulatorX.dll resource string tables 31-35 via Ghidra.
Indices 0-55. We have raw stage JSON for 33 of these 56 (indices 0-32).
Indices 33-55 (23 names) are in the name table but have no extracted stage
data in any p2k_skins dir. The user's "50 P2K skins" recollection was an
approximation of this 56-entry name space.

## Bit-exact parity gate (built this session)

Fully wired. All 35 canonical P2K skins null at **-151.9 dB** (float32
round-trip quantization floor) via `python tools/parity_null.py`.

### New tools

**`tools/render_canonical_refs.py`** (new, 100 lines):
- Reads every `*.json` in `trenchwork_clean/datasets/p2k_skins/` (35 skins).
- Pipeline: raw stage → `[b0,b1,b2,1,a1,a2]` SOS matrix → `scipy.signal.sosfilt(dry)` → `apply_agc` (16-entry table, python port of trench-core agc.rs) → `* boost`.
- Dry input: `trenchwork_clean/ref/bypassed-pinknoise.wav` (132694 samples, 44100 Hz, mono `[:,0]`).
- Writes 4 corner WAVs per skin to `Trench/ref/canonical/<name>_<corner>.wav` as float32 (subtype `FLOAT`, not PCM_16 — PCM_16 clips M100 corners whose peak exceeds 1.0 after the 4× boost).
- Writes `ref/canonical/MANIFEST.json` listing each skin with source path, boost, stage_count, rendered corners.
- Total output: 140 reference WAVs + 1 manifest.

**`tools/parity_null.py`** (rewritten, ~210 lines):
- Enumerates `ref/canonical/MANIFEST.json`, loads each skin's raw JSON from `p2k_skins/`, re-renders via the same pipeline, compares to the corresponding canonical WAV.
- Null method: `gain = dot(pred, ref) / dot(pred, pred)`, `res = ref - gain*pred`, `rel_null = 20*log10(rms(res) / rms(ref))`. Lag search disabled (`lag_search=0`) because references are rendered through the exact same pipeline — no offset to discover.
- Threshold: `FAIL_THRESHOLD_DB = -140.0`. Anything above fails with exit code 1 and lists offending skins.
- Output: one line per skin showing worst-corner null. Final `OK: N/N` or `FAIL: k/N` line.
- Uses `gc.collect()` per-skin to survive Windows low-memory conditions (system was at 95%+ RAM this session from an unrelated process; numpy was failing 1 MB allocations without the GC force).

**`ref/canonical/PROVENANCE.json`** — chain of custody pointing at canonical
source, dry input, pipeline description, sample rate, stage count, boost,
subtype. Also notes the legacy hedz*.wav files as the thing it replaces.

**`ref/canonical/MANIFEST.json`** — 35 entries, one per skin, with source
path, boost, stage_count, corner list.

### Parity gate current state

```
$ python tools/parity_null.py
dry input:  bypassed-pinknoise.wav  132694 samples
refs:       C:\Users\hooki\Trench\ref\canonical  (35 skins)
pipeline:   raw-ROM → SOS cascade → AGC → ×boost → lag+gain null
fail at:    rel_null > -140 dB

skin                      boost worst_corner   lag     gain    rel_null
 P2k_000                   4.00      M100_Q0     0   1.0000   -151.9 dB
 P2k_001                   4.00    M100_Q100     0   1.0000   -151.9 dB
 ...                       ...          ...     0   1.0000   ...
 Vocal_Oo_Ah               4.00      M0_Q100     0   1.0000   -151.9 dB

OK: 35/35 skins null at <= -140 dB
```

All 35 skins null at the float32 quantization floor. Gate passes. Exit 0.

### Null score definition

`rel_null` is the dB ratio of residual RMS to reference RMS after a
best-fit scalar gain (and optionally integer lag):
```
gain = dot(pred, ref) / dot(pred, pred)
res  = ref - gain*pred
rel_null = 20*log10(rms(res) / rms(ref))
```
Negative, lower = better. 0 dB = no nulling. -150 dB = float32 round-trip
floor. -280 dB = f64 numerical noise. -∞ = perfect. Our gate sits at -140
to leave 12 dB of headroom above the f32 round-trip floor.

## Cascade + AGC + boost chain (verified this session)

Python reference pipeline (used to render `ref/canonical/` and to verify it):

```python
def stage_co(s):
    a1 = float(s["a1"])
    r  = min(float(s["r"]), 0.999999)
    a2 = r * r
    if float(s.get("flag", 1.0)) < 0.5:
        b0, b1, b2 = 1.0, 0.0, 0.0
    else:
        b0 = 1.0 + float(s["val1"])
        b1 = a1  + float(s["val2"])
        b2 = a2  - float(s["val3"])
    return [b0, b1, b2, 1.0, a1, a2]

sos = np.array([stage_co(s) for s in stages], dtype=np.float64)
cascaded = scipy.signal.sosfilt(sos, dry)     # f64 throughout
pred = apply_agc(cascaded) * boost             # AGC per-sample, then boost
```

AGC (per-sample, matches `trench-core/src/agc.rs::agc_step` verbatim):

```python
AGC_TABLE = np.array(
    [1.0001, 1.0001, 0.996, 0.990, 0.920, 0.500, 0.200, 0.160,
     0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120, 0.120], dtype=np.float64)

def apply_agc(samples):
    out = np.empty_like(samples)
    g = 1.0
    for i in range(len(samples)):
        s = samples[i]
        idx = int(g * abs(s)) & 0xF
        ng = g * AGC_TABLE[idx]
        g = ng if ng < 1.0 else 1.0
        out[i] = s * g
    return out
```

At pink-noise levels (peak ≈ 0.23 after dry, peak ≈ 1.8 after 4× boost on
M100 corners), AGC is effectively passthrough — `int(gain*|sample|) & 0xF`
stays in [0,1] where table value is 1.0001, clamped back to 1.0 next
sample. AGC gain reduction only triggers when `|sample|*gain ≥ 2`. Our
test signal never goes there. Result: `pred` equals `cascaded * boost`
bit-exact for this dry input. Future dries that push the cascade harder
will start exercising the AGC code path, and the parity gate will
immediately catch any divergence.

Per-corner boost in canonical P2k_013 is the single top-level `boost: 4.0`.
No per-corner boost variation. All 35 skins currently use boost=4.0 (but
the renderer reads it from the skin file, so skins with different boost
values would be handled transparently).

## Commits landed this session

Branch `chore/git-hygiene`:

- `7ecfbbc session state rev 2: post-execution handoff dump`
- `42e798b audition: LP-in-morph candidate body (6-stage, no 7th global LP)` — `cartridges/auditions/lp_in_morph.json` (375 lines)
- `cfc3f70 port AGC from trenchwork_clean trench-core` — `trench-core/src/agc.rs` (89 lines, 5 unit tests pass), `trench-core/src/lib.rs` adds `pub mod agc`

Earlier in the session (committed in background): `0cb3f29`, `507bca7`,
`5a1ae92`, `8d30aae`, `adab1e9` — lean workspace, cartridge schema,
./check, parity_null v1, null target tests, research docs.

**This rev 3 block of work is UNCOMMITTED on disk:**

- `tools/render_canonical_refs.py` (new, ~120 lines)
- `tools/parity_null.py` (rewritten — removed lag search, switched to per-skin iteration with MANIFEST, added SKINS constant, added FAIL_THRESHOLD_DB, DRY+SKINS+REF_DIR constant layout)
- `ref/canonical/` (new directory, 140 wavs + MANIFEST.json + PROVENANCE.json)
- `dev/SESSION_STATE.md` (this file, rev 3)

Operator must commit these before clearing the session. Suggested commit
split:
1. `parity: canonical P2K references + 35-skin bit-exact gate` — covers
   `tools/render_canonical_refs.py`, `tools/parity_null.py`, `ref/canonical/`
2. `session state rev 3` — `dev/SESSION_STATE.md`

## Domain 1 truths reconfirmed or added this session

1. **Trench/cartridges/p2k/ compiled-v1 kernel form is direct biquad**
   (`c0..c4 = b0,b1,b2,a1,a2`), NOT the shifted minifloat-domain form used
   by `trenchwork_clean/pyruntime/encode.py`. Verified at 1e-7 against
   canonical raw for all 4 corners of P2k_013.

2. **Canonical P2K raw stage data is at `trenchwork_clean/datasets/p2k_skins/`**
   (35 skins). Byte-identical copies exist in ~15 workspace locations. All
   agree. This is the immutable source of truth for P2K filter bodies.

3. **`trenchwork_clean/cartridges/00_talking_hedz.json` has M0_Q100 ↔ M100_Q0
   stage data swapped.** Do not load it in any parity pipeline.

4. **`C:/Users/hooki/do-it/archive/data-dumps/TalkingHedz_Complete.json` is
   a broken Cheat Engine capture** with 3 independent bugs: off-by-one
   stage indexing (5 stages vs canonical 6), `val1` slot stores `1+val1`
   (b0) instead of raw val1, radius values drift 0.003-0.049 from canonical.
   The legacy `trenchwork_clean/ref/hedz*.wav` references were rendered
   from this data.

5. **Legacy `trenchwork_clean/ref/hedz*.wav` are unusable for parity.**
   Should be quarantined or deleted. Replaced by `Trench/ref/canonical/`.

6. **Parity gate passes at -151.9 dB across all 35 P2K skins.** The gate
   is a deterministic regression check for cascade math, AGC table, and
   boost application. Any future divergence in any of these produces a
   null worse than -140 dB and fails the gate.

7. **Minifloat coefficient quantization does not degrade high-Q nulls.**
   Tested round-trip encode→decode→recombine via
   `pyruntime.encode.raw_to_encoded` + `pyruntime.minifloat.pack/decode`
   against all 4 corners of the broken cartridge. M0_Q0 nulled 18 dB
   better with round-trip (because the reference WAS rendered through
   the minifloat path); M0_Q100/M100_Q0/M100_Q100 were unchanged. The
   rev 2 hypothesis that minifloat limit cycles caused the Q=100 bass
   residual is falsified. The actual cause was upstream CE capture bugs.

8. **f32 vs f64 precision does not explain the Q=100 residual.** Tested
   `sosfilt` at f32, `sosfilt` at f64, and `lfilter` (direct tf form)
   at f64. All three give identical nulls at every corner. Rules out
   state-variable precision drift.

9. **The cascade magnitude response of canonical vs reference is flat
   at 0 ± 0.1 dB across all bands at all corners** when pipeline and
   reference are aligned. The gap in rev 2 was time-domain-only, not
   spectral. Re-verified after switching to canonical refs: all nulls
   land at -151.9 dB with gain = 1.0000 and lag = 0.

## Trench-mcp tooling state

`trench-mcp` MCP server is running. Tool schemas loaded this session:
`analyze_state`, `compare_states`, `probe_frequency`, `scan_stability`,
`scan_vault`, `vault_seeds`. All operate on compiled-v1 cartridges via
`cartridge_path`. Tested `analyze_state` on
`Trench/cartridges/p2k/P2k_013.json` at morph=1, q=1 — returns terrain
(ridge/basin/saddle), shape (centroid 277 Hz / 20.4 dB), vowel-space
(talkingness 0.62, vowel bias "eh"), reactor stress 0.40 NOMINAL.
Useful for authoring-side body analysis; not currently wired into
parity or compile pipelines.

## Blockers

1. **System RAM at 95%+** this session (likely another process on the
   machine). Numpy was failing 1 MB allocations without explicit
   `gc.collect()` per-iteration in the parity loop. The loop is
   rewritten to handle this, but if the condition persists, operator
   may need to restart the shell or investigate what is holding memory.

2. **`Trench/trench-core/` still has no `engine.rs`.** `cascade.rs` and
   `agc.rs` are in place, but the full `FilterEngine` with motion
   pursuit, DC blocker, control-rate parameter dispatch, Q extrapolation
   guard, and gain ceiling clamp was not ported from
   `trenchwork_clean/trench-core/src/engine.rs` (1429 lines). Dedicated
   port session needed.

3. **No rust integration test against `ref/canonical/`.** The python
   gate works end-to-end; the rust side still has no binding. Next
   step: small rust binary at `trench-core/tests/canonical_parity.rs`
   (or similar) that loads each `Trench/cartridges/p2k/*.json`, runs
   the frozen `Cascade` + ported `agc_step` + boost, writes corner
   WAVs to `ref/rust_render/`, and asserts bit-exact null (within f32
   tolerance) vs `ref/canonical/`. This is the bridge from python
   reference to rust production.

4. **Hardware benchmark captures remain orphaned.** 33 files at
   `trenchwork_clean/captures/benchmark/P2K_*.wav` (8 s mono 44100 Hz,
   plausibly real hardware recordings via sine sweeps) still have no
   located dry input. The short pink-noise dry
   (`bypassed-pinknoise.wav`, 3 s) does not match. Options to unlock
   them: (a) find a matching sweep input file; (b) characterize the
   sweep parameters from one capture and generate a matching stimulus;
   (c) sweep-deconvolve the captures to impulse responses and null
   IR-vs-IR. Not blocking the current python+rust parity chain.

5. **Branch is still `chore/git-hygiene`**, not `codex/v1`. This
   session's commits are accumulating on a side branch. Merge/rebase
   strategy still not decided.

6. **Hardware parity step in `./check` is report-only.** Previously
   printed numbers against broken hedz*.wav; now that parity_null.py
   is bit-exact against canonical refs, `./check` should be updated
   to FAIL on non-zero exit from parity_null. One-line change in
   `./check`; not yet applied.

## Next steps (ordered)

1. **Commit rev 3 work.** Split into two commits as noted above.
2. **Update `./check`** to propagate the parity_null.py exit code as a
   hard gate. One-line change: remove the `|| true` or equivalent
   soft-skip if present; otherwise already correct.
3. **Port `FilterEngine` + DC blocker** from
   `trenchwork_clean/trench-core/src/engine.rs` (1429 lines) into
   `Trench/trench-core/src/`. Keep the Trench `Cascade` frozen per
   doctrine; the port is a new `engine` module built on top of the
   existing `Cascade` + the already-ported `agc`. Dedicated session.
4. **Write rust `canonical_parity` integration test.** Loads each
   `cartridges/p2k/*.json`, processes the dry pink noise via the rust
   cascade + AGC + boost, writes corner WAVs to `ref/rust_render/`,
   asserts bit-exact null vs `ref/canonical/` (threshold ≈ -140 dB to
   match the python gate, or tighter if rust does all-f64 internally).
   This binds the python reference chain to the rust shipping path.
5. **Delete or quarantine `trenchwork_clean/ref/hedz*.wav`.** The files
   are rendered from the broken CE capture and have no further purpose.
   Keep one copy tagged `_BROKEN_CE_CAPTURE_QUARANTINE` if the operator
   wants it preserved as evidence. Update any docs that reference these
   wavs.
6. **Locate or reconstruct dry input for benchmark captures.** Next-tier
   goal — unlocks hardware-level parity testing for all 33 P2K presets
   against real hardware recordings. Not blocking v1 shipping.
7. **Author Speaker Knockerz / Cul-De-Sac / Aluminum Siding / Small Talk
   Ah-Ee candidates** using the LP-in-morph pattern proven in
   `cartridges/auditions/lp_in_morph.json`. Each candidate validated
   against `cartridge.schema.json` and its `BODIES.md` "how to know it's
   wrong" rubric.
8. **Decide branch strategy** for `chore/git-hygiene` — merge to
   `codex/v1`, rebase, or open PR. Operator decision.
9. **Fill name coverage for indices 33-55 in `p2k_filter_names.json`.**
   23 named P2K filters are in the EmulatorX.dll resource string table
   but have no extracted raw-ROM stage JSON in `p2k_skins/`. If those
   filters exist in the ROM at addresses not yet parsed, add them to
   the canonical set. If they don't exist in the ROM (string-table-only
   slots), document that explicitly.

## Operator handoff cue

When the operator resumes, the python parity chain is complete and
passing. The cleanest next action is **commit rev 3 work** (step 1),
then start the **rust FilterEngine port** (step 3) which is the bridge
to rust parity. The python gate gives an immediate regression check
while porting — any rust test that produces WAVs differing from
`ref/canonical/` by more than the f32 quantization floor is evidence
of a cascade/AGC/boost divergence that needs immediate attention.

The rev 2 "Q=100 bass residual" and "M0_Q0 2× amplitude" mysteries are
both fully closed and should not consume any further analysis time.
Upstream Cheat Engine capture bugs (off-by-one stages, val1-as-b0,
radius drift) were the root cause. The canonical 6-stage raw ROM
extraction has been the correct source all along; the parity pipeline
now uses it exclusively.
