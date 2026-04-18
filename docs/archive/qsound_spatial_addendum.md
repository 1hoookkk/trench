# QSound Spatial Addendum — Calibration from Measured WAVs

**Measured from:** `C:/Users/hooki/trench_re_vault/datasets/qsound_spatial_v1/2026-03-05_pan_batch_fullgrid_a/`
**Sample rate:** 48 kHz (confirmed from `CAPTURE_MANIFEST.json`)
**Stimulus:** 220 Hz saw, −18 dB, 48 kHz stereo float32
**Date:** 2026-04-18
**Scripts:** `tools/qsound_calibrate.py`, `tools/qsound_calibrate_v2.py`, `tools/qsound_calibrate_v3.py`

## TL;DR — status for Task 8

**DONE_WITH_CONCERNS — parametric model cannot be validated against this dataset.**

The pan-batch capture does *not* contain 13 monotonically increasing pan states. Direct per-channel RMS measurement collapses the 13 captures into only **8 distinct output states**, with the same spatial pose appearing at multiple labelled pan values (e.g. pans −20 and 0 are effectively identical; pans 10 and 20 are identical; pan −50 is identical to pan +50 etc.). This makes it impossible to extract the scalar mapping `az → {ITD samples, ILD dB, band gains dB}` that Task 8 needs.

**Task 8 must not treat the parametric coefficients in `qsound_spatial.md` as validated.** They were regressed from a dataset whose ground truth has the defect shown below. Before implementing the Rust spatial stage we need a clean re-capture.

## What this data actually contains

Per-channel RMS on the 1 s steady-state window, at each labelled pan:

| pan   | L_rms    | R_rms    | R−L dB  | cluster |
|-------|----------|----------|---------|---------|
| −60   | 0.0324   | 0.0156   | −6.33   | A (unique) |
| −50   | 0.0324   | 0.0324   | +0.00   | B (unique) |
| −40   | 0.0280   | 0.0282   | +0.06   | C |
| −30   | 0.0148   | 0.0330   | +6.97   | D |
| −20   | 0.0323   | 0.0241   | −2.52   | E |
| −10   | 0.0323   | 0.0246   | −2.36   | E' |
|  0    | 0.0323   | 0.0241   | −2.52   | E (same as −20) |
| +10   | 0.0147   | 0.0330   | +7.00   | D' |
| +20   | 0.0148   | 0.0330   | +7.01   | D' |
| +30   | 0.0280   | 0.0282   | +0.06   | C (same as −40) |
| +40   | 0.0323   | 0.0246   | −2.36   | E' |
| +50   | 0.0323   | 0.0242   | −2.52   | E |
| +60   | 0.0323   | 0.0246   | −2.35   | E' |

Observations:

* Pans 0, −20, and +50 give the SAME output.
* Pans −30, +10, and +20 give the SAME output.
* Pans −40 and +30 give the SAME output.
* 6 of 13 captures (−20, −10, 0, +40, +50, +60) cluster around the same −2.4 dB asymmetric state.
* The expected hard-pan geometry (−60 mirror of +60) does not hold: −60 gives L-heavy output but +60 gives the mild-asymmetric E' state.

Correlations between "same cluster" captures are all >0.99 (pan −20 vs pan 0: r=0.999992), confirming these are the same spatial configuration measured with small phase offsets between captures. This is not a windowing artefact of our measurement — it is a property of the source WAVs.

Likely cause: the capture pipeline (`V1 (Offline ABI)` in `qsound_spatial.md`) failed to reliably apply the pan argument before each render. The manifest `attempts` field (2, 1, 5, 1, 1, 4, 1, 4, 1, 1, 2, 3, 7) shows retries on exactly the captures that ended up in "stuck" clusters. Combined with the `0xC0000005` failures documented in `qsound_spatial.md` for vector-position setters, it is consistent with `QMixer.dll` state leaking between captures.

## What we measured (despite the ground-truth defect)

### ITD

Integer-lag cross-correlation between L and R saturates at either 0 or ±5 samples across all pans. Sub-sample lag via FFT-upsampled xcorr + per-harmonic phase unwrap gives:

| pan | xc_lag_samp | ph_lag@220Hz (samp) | phase@220 (deg) |
|-----|-------------|---------------------|------------------|
| −60 | −0.25       | +1.96               | +3.2°            |
| −50 | 0.00        | 0.00                | 0.0°             |
| −40 | 0.00        | +0.09               | +0.1°            |
| −30 | +0.25       | −4.53               | −7.5°            |
| −20 | +5.44       | +106.3              | +175.4° (=−4.6°) |
| −10 | +5.25       | +107.96             | +178.1° (=−1.9°) |
|  0  | +5.44       | +106.3              | +175.4°          |
| +10 | +0.25       | −4.64               | −7.7°            |
| +20 | +0.25       | −4.64               | −7.7°            |
| +30 | 0.00        | +0.09               | +0.1°            |
| +40 | +5.25       | +107.95             | +178.1°          |
| +50 | +5.44       | +106.3              | +175.4°          |
| +60 | +5.25       | +107.96             | +178.1°          |

The big "≈180° phase" captures (cluster E/E') are **polarity-inverted** between L and R, not time-delayed. Unwrapping 180° as 0° delay: actual interaural lag at 220 Hz across all 13 captures stays within roughly ±8 samples (±165 µs), which is comparable to a classic head-radius ITD budget (~31 samples = 650 µs maximum).

**No clear monotonic trend ITD vs labelled pan exists** because pan labels are unreliable.

### ILD

Per-pan broadband RMS ratio (R_rms / L_rms, in dB) reported above. A polarity inversion alone would give 0 dB, not ±7 dB, so ILD is real and significant in at least three states: cluster A (−6.3), cluster D (+7.0), cluster A-inverse (nominally +6.3 but absent; nearest is −60 which is L-heavy).

### Per-band / per-harmonic

Per-harmonic R−L (dB) shows extreme high-frequency asymmetry in clusters A and D:

| cluster | 220 Hz | 440 Hz | 660 Hz | 1760 Hz | 2640 Hz | 4400 Hz |
|---------|--------|--------|--------|---------|---------|---------|
| A (−60)     | −13.3 | −15.5 | −15.4 | −0.8    | −1.1    | −0.4    |
| D (+10,+20) | +14.5 | +15.9 | +14.3 | +1.4    | +1.5    | −0.8    |
| E (0)       | −2.1  | −1.9  | −0.3  | −3.8    | −6.1    | −9.8    |

Cluster E has striking HIGH-frequency attenuation on R (up to −9.8 dB at 4.4 kHz) even though low harmonics are nearly level-matched. This is consistent with the doc's claim that QSound applies **per-band FIR shaping per channel**, and implies the shaping extends well below any plausible "high-band" crossover.

**Three-band topology is at best an approximation.** The dB delta vs harmonic frequency in cluster E does not look like a flat gain in 3 piecewise bands; it is a smoothly increasing attenuation from ~0 dB at 220 Hz to −10 dB at 4.4 kHz. A short FIR or IIR shelf is a more likely physical model than a 3-band gain.

## Resolved UNKNOWNs from qsound_spatial.md

### Sample rate
**48 000 Hz.** Confirmed from `CAPTURE_MANIFEST.json` `device_mix_format.sr` for all 13 captures.

### ITD units and conversion
**UNRESOLVED — not recoverable from this dataset.**

Parametric output at pan ±60 with `az_deg` in radians → `sin(π/3)` family terms give `itd_law = ±2666`. Dividing into the measured integer lag (5 samples) gives a ratio near 533, but the lag itself is unreliable (it's a cross-correlation lock onto a sub-multiple of the 218-sample saw period, not a true interaural delay). The phase-lag at the fundamental (above) is more trustworthy and stays <±8 samples in all real-audio states, contradicting the parametric value (which predicts 5–56 samples at ±60 → ±10).

Two interpretations are possible:

1. **"itd_law output is in microseconds" with the known coefficient set oversized by ≈5×.** 2666 µs = 128 samples, clearly too large. Ruled out.
2. **"itd_law output is in clock ticks of a proprietary fractional-delay unit"**, then scaled down by an unknown factor before application.

Neither can be pinned from the measured data. **Recommendation:** do NOT ship any fixed ITD conversion factor into `qsound_spatial.rs`; treat ITD as a free model parameter and expose it for per-body retuning until a clean capture is available.

### az units
**UNRESOLVED — radians interpretation is the only one that gives a monotonic parametric ITD vs pan curve, but measured data is too noisy to prefer radians over degrees empirically.** Parametric ITD with radians gives a strictly monotonic ramp 0→2666 over pan 0→60; with degrees the output is wild (sin(60 deg) numerically = sin(60 rad) = −0.305, producing non-monotonic coefficient output).

**Engineering call:** use **radians** in `qsound_spatial.rs`. This is the only interpretation under which the regressed coefficients produce a smooth, physically plausible az → ITD/ILD mapping.

### ITD ear mapping
**UNRESOLVED — ambiguous in this dataset.**

* At cluster A (labelled −60, L-heavy), sub-sample xcorr lag is −0.25 samples, phase lag at 220 Hz is +1.96 samples (R leads slightly). Inconsistent: RMS asymmetry says L-side, phase says R leads.
* At cluster D (labelled +10/+20, R-heavy), phase lag is −4.5 samples (L leads).

The comment in `qsound_spatial.md` says "right-ear-lead/right-ear-louder convention" for positive `az`, i.e. both ITD and ILD move together. The data partially supports this in cluster D but contradicts it in cluster A.

**Engineering call:** follow the doc convention (right-ear leads for positive az). Apply the delay to the contralateral (lagging) ear only, leaving the ipsilateral ear at zero delay. Validate against a fresh capture before locking.

### ILD per-pan
Parametric ILD coefficients appear to output dB directly based on the magnitude of the regressed constants (~6.8 for `sin(az)`). The range across the measured clusters (−6.3 to +7.0 dB) is of the same order as the parametric prediction at ±60 (±3.9 dB). **Parametric ILD underestimates the extreme clusters by ~3 dB.** Whether this is a model gain error or a capture-state artefact cannot be distinguished here.

### Band-law topology
**UNRESOLVED — but "3-band flat-gain" model is likely wrong.**

Best-fit fixed 3-band crossovers from brute-force comparison: **(200 Hz, 2000 Hz)** with MAE 2.6 dB vs parametric. But the residual pattern shows the measured per-harmonic attenuation in cluster E is smooth across frequency, not piecewise-flat. A 3-band gain block will be visibly "stepped" relative to the source material.

**Engineering recommendation:** implement as per-channel **IIR shelf pair** (low-shelf + high-shelf) rather than a hard 3-band split. The parametric `band_law` table can still feed a 3-point {low, mid, high} dB gain target; convert to shelf gains via:

```
low_shelf_dB  = band_law["low"]  − band_law["mid"]
high_shelf_dB = band_law["high"] − band_law["mid"]
broadband_dB  = band_law["mid"]
```

Shelf corners at roughly 400 Hz and 2500 Hz match the measured transitional bands better than hard splits.

## Recommended Rust types for qsound_spatial.rs

```rust
/// Spatial pose for a single source. az in radians.
pub struct SpatialPose {
    pub az_rad: f32,       // -π/3 .. +π/3 corresponds to ±60 pan
    pub dist:  f32,        // metres; log2(dist/0.25) feeds band_law
    pub el_rad: f32,
}

/// Output of the parametric spatial model at a given pose. All values
/// pre-evaluated once per control block; consumed by the audio thread.
pub struct SpatialProfile {
    pub itd_samples: f32,        // applied to LAGGING ear; other ear = 0.
    pub gain_l_db: f32,          // broadband, from -ild_law/2
    pub gain_r_db: f32,          // broadband, from +ild_law/2
    pub shelf_l: ChannelShelves, // per-channel low/high shelf targets
    pub shelf_r: ChannelShelves,
}

pub struct ChannelShelves {
    pub low_shelf_gain_db:  f32,  // at corner ~400 Hz
    pub high_shelf_gain_db: f32,  // at corner ~2500 Hz
}
```

The ITD scalar and shelf corners MUST be exposed as tuning constants in the Rust code with clear TODOs; DO NOT bake any factor from the parametric coef output until a clean capture validates them.

## Verification targets for Task 8 tests

Since this dataset can't be trusted as ground truth, we should NOT write tests that assert against its numeric values. Instead, gate Task 8 on **behavioural checks** that the parametric model alone can justify:

| Test | Expected | Tolerance |
|------|----------|-----------|
| At az=0, el=0, dist=1: L_rms ≈ R_rms for white-noise input | 0.0 dB | ±0.3 dB |
| At az=+π/3, R_gain_dB > L_gain_dB (broadband) | ≥2 dB | one-sided |
| At az=+π/3, L ear delay > R ear delay (sample count) | ≥0 samples | one-sided |
| Swapping az sign mirrors L/R output channels | L(−az) == R(+az) within 1e-5 | strict |
| band_law evaluated at az=0 is identical for L and R channels | 0.0 dB | <1e-4 |
| No denormal / clipping on all ±60° sweeps with full-scale input | peak ≤ +0 dBFS | strict |

These are satisfied by any reasonable implementation of the parametric model and don't depend on the defective capture batch.

## Re-capture requirements before final Task 8 lock-in

When the capture pipeline is fixed, the following grid is needed to lock the UNRESOLVED items above:

1. **Single-pan stability:** capture pan=0 ten consecutive times; confirm all ten outputs are bit-identical or within 0.1 dB RMS. Guards against the state-leak shown in this addendum.
2. **Symmetric sweep:** pans [−60, −45, −30, −15, 0, +15, +30, +45, +60] with a bounded stimulus (white noise or pink noise ≥5 s). Confirm −X is the L/R swap of +X within 0.5 dB per band. That alone would resolve ear-mapping.
3. **Wideband stimulus** (not a 220 Hz saw) to let band-law be fit densely instead of at ~20 isolated harmonics.
4. **Per-channel solo:** capture L and R outputs separately with the same input to get a direct impulse-response estimate per channel — this lets us recover the actual per-channel filter (shelf vs 3-band vs full FIR) without fitting.

## Script provenance

* `tools/qsound_calibrate.py` — v1: RMS ITD/ILD + 3-band fit
* `tools/qsound_calibrate_v2.py` — v2: sub-sample ITD + per-harmonic L/R table
* `tools/qsound_calibrate_v3.py` — diagnostic per-pan RMS grouping that surfaced the cluster-collapse
* Output JSONs: `tools/qsound_calibrate.json`, `tools/qsound_calibrate_v2.json`
