# DSP Spec — FIELD (Trench)

> Moved from CLAUDE.md to keep Claude prompts lightweight. Reference this file only when working on DSP details.

## Identity

FIELD is a 14-pole Z-Plane filter effect. It generates the "morphing formant" sound of 1990s E-mu hardware (Morpheus, Audity 2000) using pure math. No samples. No ROM dumps. Coefficients calculated in real-time from pole positions.

Target: Underground producers. Gatekept sound design. Industrial/digital aesthetic.

---

## The Engine

### Data Structures
```cpp
struct Pole {
    double freq;      // Hz (20 - 20000)
    double radius;    // 0.0 - 0.999
    bool active;      // false = bypass this stage
};

struct FilterFrame {
    Pole poles[7];
    Topology topology;  // PARALLEL or CASCADE
};

struct Biquad {
    double b0, b1, b2;  // Numerator (zeros)
    double a1, a2;      // Denominator (poles)
};
```

### Structure
```
7 Biquad Stages × 2 Poles = 14 Poles Total

Each biquad: y[n] = b0·x[n] + b1·x[n-1] + b2·x[n-2] - a1·y[n-1] - a2·y[n-2]
```

### Topology Modes
```
PARALLEL: Sum outputs of all stages → Vowels, formants, EQ
CASCADE:  Chain stages in series   → Phasers, flangers, LP/HP
```

### Precision
- Coefficient calculation: `double` (64-bit)
- Audio processing: `float` (32-bit)
- State variables: `double` (prevents denormals)

---

## The Math

### Core Formula: Pole to Coefficients
```cpp
// Given:
//   freq = center frequency (Hz)
//   bw   = bandwidth (Hz) — smaller = more resonant
//   sr   = sample rate (Hz)

// Step 1: Calculate pole position
double θ = 2.0 * PI * freq / sr;        // Angle (frequency)
double R = exp(-PI * bw / sr);           // Radius (resonance)

// Step 2: Denominator (same for all types)
double a1 = -2.0 * R * cos(θ);
double a2 = R * R;
```

### The "Ring" Formula
```
Bandwidth 60-100 Hz → R ≈ 0.99  → Q ≈ 50   → Gentle
Bandwidth 15-25 Hz  → R ≈ 0.998 → Q ≈ 300  → E-MU RING ✓
```

### Coefficient Tables by Filter Type

**RESONATOR (Vowels, Formants)** — Use with PARALLEL topology
```cpp
double gain = 1.0 - R;
b0 = gain;
b1 = 0.0;
b2 = -gain;
a1 = -2.0 * R * cos(θ);
a2 = R * R;
```

**ALLPASS (Phasers, Flangers)** — Use with CASCADE topology
```cpp
b0 = R * R;
b1 = -2.0 * R * cos(θ);
b2 = 1.0;
a1 = -2.0 * R * cos(θ);  // Same as b1
a2 = R * R;               // Same as b0
```

**LOWPASS** — Use with CASCADE topology
```cpp
double K = tan(PI * freq / sr);
double norm = 1.0 / (1.0 + K/Q + K*K);
b0 = K * K * norm;
b1 = 2.0 * b0;
b2 = b0;
a1 = 2.0 * (K*K - 1.0) * norm;
a2 = (1.0 - K/Q + K*K) * norm;
```

**BYPASS (Unused stage)**
```cpp
b0 = 1.0;
b1 = 0.0;
b2 = 0.0;
a1 = 0.0;
a2 = 0.0;
```

---

## Morphing

**NEVER interpolate coefficients (a1, a2, b0, b1, b2) directly.**

Interpolating raw coefficients crosses outside the stability triangle and causes filter explosions.

**ALWAYS interpolate pole parameters, then regenerate coefficients.**

---

### Interpolation Methods

**RADIUS: Geometric (log) interpolation**
```cpp
// Perceptually linear resonance change
// Decay rate is linear, ring time feels even
double R = pow(R_A, 1.0 - t) * pow(R_B, t);
```

**FREQUENCY: Logarithmic (pitch) interpolation**
```cpp
// Perceptually linear pitch sweep
// Octaves feel evenly spaced
double freq = freq_A * pow(freq_B / freq_A, t);

// Clamp to valid range
freq = clamp(freq, 20.0, sampleRate * 0.45);
```

**ANGLE: Never wrap**
```cpp
// Treat frequency as line segment [0, fs/2], not circle
// 18kHz → 200Hz sweeps DOWN through the spectrum
double θ = 2.0 * PI * freq / sampleRate;
```

---

### Pole Assignment (when stage counts differ)

When Frame A has 4 poles and Frame B has 7, you must MATCH poles to prevent "teleporting."

**Cost function (perceptual distance):**
```cpp
double cost(Pole a, Pole b) {
    double d_freq = abs(log(a.freq) - log(b.freq));
    double d_radius = abs(log(a.radius) - log(b.radius));
    return d_freq + d_radius;
}
```

**Assignment:** Brute force all 7! = 5040 permutations, pick minimum total cost.
Do this ONCE when setting Frame A/B, not every 32 samples.

**Unmatched poles → Birth/Death:**
- Birth (inactive→active): Morph FROM bypass coefficients TO designed coefficients
- Death (active→inactive): Morph FROM designed coefficients TO bypass

```cpp
// Birth: pin frequency to destination, fade in
if (birth) {
    // Frequency stays at target (no "sweep in from nowhere")
    freq = frameB.poles[i].freq;
    // Radius fades from 0 (bypass) to target
    R = pow(frameB.poles[i].radius, t);  // 0→target
}

// Death: pin frequency to source, fade out  
if (death) {
    freq = frameA.poles[i].freq;
    R = pow(frameA.poles[i].radius, 1.0 - t);  // target→0
}
```

---

### Topology Changes (PARALLEL ↔ CASCADE)

**Cannot morph structure. Must crossfade outputs.**

```cpp
if (frameA.topology != frameB.topology) {
    // Run BOTH filters during transition
    double yParallel = processParallel(x);
    double yCascade = processCascade(x);
    
    // Crossfade outputs
    return (1.0 - t) * yParallel + t * yCascade;
}
```

This is always stable (sum of two stable outputs).

---

### Complete Morph Update (every 32 samples)

```cpp
void updateMorph(float t) {  // t = 0.0 to 1.0
    for (int i = 0; i < 7; ++i) {
        // Skip if both inactive
        if (!frameA.poles[i].active && !frameB.poles[i].active) {
            setBypass(i);
            continue;
        }
        
        // Birth: inactive → active
        if (!frameA.poles[i].active && frameB.poles[i].active) {
            double freq = frameB.poles[i].freq;
            double R = pow(frameB.poles[i].radius, t);
            regenerateCoeffs(i, freq, R);
            continue;
        }
        
        // Death: active → inactive
        if (frameA.poles[i].active && !frameB.poles[i].active) {
            double freq = frameA.poles[i].freq;
            double R = pow(frameA.poles[i].radius, 1.0 - t);
            regenerateCoeffs(i, freq, R);
            continue;
        }
        
        // Normal morph: both active
        // Log interpolation for frequency
        double freq = frameA.poles[i].freq * 
                      pow(frameB.poles[i].freq / frameA.poles[i].freq, t);
        
        // Geometric interpolation for radius
        double R = pow(frameA.poles[i].radius, 1.0 - t) * 
                   pow(frameB.poles[i].radius, t);
        
        regenerateCoeffs(i, freq, R);
    }
}

void regenerateCoeffs(int i, double freq, double R) {
    double θ = 2.0 * PI * freq / sampleRate;
    coeffs[i].a1 = -2.0 * R * cos(θ);
    coeffs[i].a2 = R * R;
    
    if (topology == PARALLEL) {
        double gain = 1.0 - R;
        coeffs[i].b0 = gain;
        coeffs[i].b1 = 0.0;
        coeffs[i].b2 = -gain;
    } else {
        coeffs[i].b0 = R * R;
        coeffs[i].b1 = coeffs[i].a1;
        coeffs[i].b2 = 1.0;
    }
}

void setBypass(int i) {
    coeffs[i] = {1.0, 0.0, 0.0, 0.0, 0.0};
}
```

---

## Saturation

The E-mu "scream" comes from saturating the FEEDBACK path, not the output.

**Dattorro Soft Clipper:**
```cpp
inline double saturate(double x) {
    if (x > 1.0) return 0.666666;
    if (x < -1.0) return -0.666666;
    return x - (x * x * x) / 3.0;  // x - x³/3
}
```

**Application (inside biquad):**
```cpp
double y = b0*x + b1*x1 + b2*x2 - a1*y1 - a2*y2;
y1 = saturate(y);  // Saturated value goes into feedback state
y2 = y1_prev;
```

---

## Gain Staging (E-mu Architecture)

**Why narrow resonators sound quiet:** With R≈0.998 (18Hz bandwidth), each formant passes only a tiny fraction of white-noise power. The peaks look right on the analyzer, but the output is inaudible without **bandwidth loudness compensation**.

**E-mu's solution:** Treat gain as a first-class parameter. Interpolate encoded coefficients, then apply:
1. Decoded filter gain
2. Aggregate gain factor
3. Soft clipping

---

### Complete Signal Flow (Parallel Formant Bank)

```
Input
  │
  ├─► (A) PreGain / Drive
  │       x0 = x * preGain * drive
  │       preGain = 1.0 (don't kneecap input)
  │
  ▼
┌─────────────────────────────────────────────────┐
│  (B) Per-Formant Resonators (parallel)          │
│                                                 │
│  Each stage i fed by SAME x0:                   │
│    θ = 2π·freq/sr                               │
│    a1 = -2R·cos(θ), a2 = R²                     │
│    b0 = g·(1-R), b1 = 0, b2 = -g·(1-R)          │
│    (zeros at DC and Nyquist)                    │
│                                                 │
│  Saturation INSIDE state path:                  │
│    y1 = saturate(y)  // before storing          │
└─────────────────────────────────────────────────┘
  │
  ├─► (C) Per-Formant Weights
  │       y_weighted[i] = y[i] * weight[i]
  │       (This is where vowel character lives)
  │
  ├─► (D) Power Normalization (not divide-by-N)
  │       y_sum = Σ(w[i] * y[i])
  │       y_norm = y_sum / (sqrt(Σ w[i]²) + ε)
  │
  ├─► (E) Bandwidth Loudness Compensation
  │       B_w = -ln(R) * sr / π  (≈30Hz for R=0.998 @ 48k)
  │       B_ref = 300 Hz (vocal reference)
  │       g_bw = sqrt(B_ref / B_w)  (≈3.16 = +10dB)
  │
  ├─► (F) PostGain + Soft Clip
  │       y = softclip(y_norm * g_bw * postGain)
  │       postGain ≈ 4.0 (+12dB)
  │
  ▼
Output (trim 0.5-1.0)
```

---

### Bandwidth Loudness Compensation Formula

```cpp
// The missing piece that makes narrow resonators audible
double computeBandwidthGain(double R, double sampleRate, double B_ref = 300.0) {
    // Bandwidth from radius: B_w = -ln(R) * sr / π
    double B_w = -std::log(R) * sampleRate / PI;

    // Constant-loudness compensation
    return std::sqrt(B_ref / B_w);
}

// For R=0.998 @ 48kHz:
// B_w ≈ 30 Hz
// g_bw = sqrt(300/30) = 3.16 (+10 dB)
```

---

### Resonator Numerator (Zeros at DC/Nyquist)

```cpp
// Zeros at z=+1 (DC) and z=-1 (Nyquist)
// B(z) = g(1-R)(1 - z⁻²)
b0 = g * (1.0 - R);
b1 = 0.0;
b2 = -g * (1.0 - R);

// g is the per-stage loudness control
// NOT "normalization" - it's an explicit gain parameter
```

---

### Per-Formant Weights

```cpp
// Vowel character comes from RELATIVE formant amplitudes
// Not just frequencies
struct FormantPreset {
    double freq;
    double bandwidth;
    double weight;  // Linear amplitude (0.5 - 2.0 typical)
};

// Example: "Ah" might weight F1 higher than F4
FormantPreset ah_formants[4] = {
    {730,  18, 1.0},   // F1 - strongest
    {1090, 18, 0.8},   // F2
    {2440, 18, 0.5},   // F3
    {3500, 18, 0.3},   // F4 - weakest
};
```

---

### Complete Gain Chain (Concrete Numbers)

For 4 formants @ R=0.998, sr=48kHz:

| Stage | Gain Factor | dB | Cumulative |
|-------|-------------|-----|------------|
| PreGain | 1.0 | 0 | 0 dB |
| Per-formant g | 1.0 | 0 | 0 dB |
| Power norm | ÷2 | -6 | -6 dB |
| BW comp | 3.16 | +10 | +4 dB |
| PostGain | 4.0 | +12 | +16 dB |
| Soft clip | ~0.8 | -2 | +14 dB |
| Trim | 0.7 | -3 | **+11 dB** |

This +11dB total lift is what makes the vowel audible when you're only passing ~120Hz total bandwidth of a 24kHz noise band.

---

## Filter Presets

### Vowels (PARALLEL topology, tight bandwidth)

| Vowel | F1 | F2 | F3 | F4 | BW |
|-------|-----|------|------|------|-----|
| Ah | 730 | 1090 | 2440 | 3500 | 18 |
| Ee | 270 | 2290 | 3010 | 3700 | 16 |
| Oo | 300 | 870 | 2240 | 3500 | 15 |
| Eh | 530 | 1840 | 2480 | 3500 | 18 |
| Uh | 640 | 1190 | 2390 | 3500 | 20 |

Stages 5-7: Set to BYPASS or wide "body" resonances (150Hz, 200Hz, 250Hz with BW=400).

### Phasers (CASCADE topology)

| Type | Freq Spread | Stages | R |
|------|-------------|--------|-------|
| Deep | 200-4000 Hz log | 7 | 0.95 |
| Extreme | 100-8000 Hz log | 7 | 0.98 |
| Subtle | 500-3000 Hz log | 4 | 0.90 |

```cpp
// Logarithmic frequency spread
for (int i = 0; i < stages; ++i) {
    float t = (float)i / (stages - 1);
    float freq = startFreq * pow(endFreq / startFreq, t);
    poles[i] = { freq, radius };
}
```

### Acid (CASCADE topology, tight single peak)
```cpp
// 4 poles stacked at cutoff for 24dB/oct slope
float cutoff = 40.0 * pow(400.0, morph);  // 40Hz - 16kHz
float R = 0.98;  // High resonance
for (int i = 0; i < 4; ++i) {
    poles[i] = { cutoff, R };
}
// Stages 5-7: BYPASS
```

---

## User Parameters

| Parameter | Range | Maps To |
|-----------|-------|---------|
| MORPH | 0-100% | Interpolation between Frame A and Frame B |
| Q | 0-100% | Global radius scaling (0.7 → 0.999) |
| DRIVE | 0-100% | Input gain into saturators |

### Q Implementation
```cpp
void applyQ(float qKnob) {  // 0.0 to 1.0
    float minR = 0.7;
    float maxR = 0.999;
    for (auto& pole : poles) {
        pole.radius = minR + qKnob * (maxR - minR);
    }
}
```

---

## Signal Flow

```
Input
  │
  ├─► Drive (user parameter)
  │
  ├─► Input Gain (calculated from pole radii)
  │
  ▼
┌─────────────────────────────────────┐
│         7 × Biquad Stages           │
│                                     │
│  PARALLEL: sum all outputs          │
│  CASCADE: chain in series           │
│                                     │
│  Saturation in feedback path        │
└─────────────────────────────────────┘
  │
  ▼
Output
```

---

## UI: Spectrum Excavator

Visual: Vertical waterfall spectrogram. SDR/hacker aesthetic.
- Dark background, toxic green signal
- Filter poles appear as bright vertical "scars"
- CRT scanlines, chromatic aberration, pixel grit
- NO KNOBS. NO SLIDERS. NO MENUS.

**Interactions:**
```
DRAG LEFT/RIGHT  = Morph (A → B)
SCROLL WHEEL     = Q (resonance)
DOUBLE-CLICK     = Resonance bloom (Q spike)
```

**Display shows:**
- Input spectrum (dim)
- Filter response (bright scars where poles are)
- Output result (where they intersect)

---

## Implementation Order

1. ZPlaneCore — 7 biquads with saturation, both topologies
2. ZPlaneDesigner — pole-to-coefficient math
3. One vowel preset — verify "Ah" sounds correct through white noise
4. Morph — interpolate between Ah and Ee, verify stability
5. Phaser preset — verify allpass coefficients work
6. Q parameter — global radius scaling
7. Basic UI — waterfall spectrogram
8. Interactions — drag for morph, scroll for Q

---

## DO NOT

- Interpolate coefficients directly (causes explosions)
- Use standard tanh (use x - x³/3)
- Clamp resonance to "safe" values (let it ring)
- Build GUI before DSP is verified
- Hardcode coefficient tables (calculate from freq/bw)
- Add safety limiters (artifacts are features)

## DO

- Calculate everything from frequency and bandwidth
- Saturate the feedback path
- Test with white noise — you should hear vowels
- Verify morph stability across full 0-100% range
- Use 64-bit doubles for coefficient math
- Make it scream when Q is maxed