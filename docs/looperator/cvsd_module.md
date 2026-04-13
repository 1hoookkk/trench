# Source
- `C:\Users\hooki\trench_re_vault\extract\cvsd\2026-03-05\cvsd_behavioral_spec.v1.yaml`

# What this data is
This captures the operational behavioral specifications for the Continuously Variable Slope Delta (CVSD) modulation block. Instead of ripping raw instruction lines, this isolates the actual deterministic acoustic response metrics (e.g. tracking slew-rates via syllabic charge) when the codec operates under test signals such as steps, impulses, and white noise. 

# Why it's here
In the CLEAN DSP codebase, the CVSD block must provide an authentic degradation path echoing early E-mu lo-fi grit. Hard-coding arbitrary saturation isn't historically accurate. Relying on these exact decay constants and idle tones ensures the runtime implementation breathes correctly behind the main filter logic inside the Looperator UI.

# Sanitized content

## Oracle Characteristics
- **Reference Oracle Tool:** SoX (Sound eXchange) `v14.4.2`
- **Canonical Sample Rate:** 16,000 Hz.
- **Determinism Status:** Certified `True`. 

## Defined Acoustic / Behavioral Parameters

| Parameter | Measured Value | Semantic Role |
|---|---|---|
| `syllabic_charge_tau_ms` | 39.499 ms | Envelope rise speed defining tracking adaptability to hard attacks. |
| `syllabic_discharge_tau_ms` | 16.546 ms | Envelope drop speed guiding compressor-like release after transients. |
| `reconstruction_integrator_tau_ms` | 4.360 ms | Leaky integrator characteristic defining the lowpass curve of the delta-encoded steps. |
| `idle_oscillation_freq_hz` | 2596.0 Hz | Background carrier idle pattern pitch. |
| `idle_oscillation_level_dbfs` | -91.29 dBFS | Quantization hum floor. |
| `transient_overshoot_db` | 2.47 dB | Clipping/ceiling bump during step-input tracking lag. |
| `transient_recovery_ms` | 2.3125 ms | Time-to-stabilize post transient snap. |
| `transient_smear_ms` | 0.0 ms | Strict zero pre/post-ringing outside of tracking error. |

**Code Representation (Behavioral YAML):**
```yaml
measured:
  syllabic_charge_tau_ms: 39.499355029870756
  syllabic_discharge_tau_ms: 16.54656555730313
  reconstruction_integrator_tau_ms: 4.360179593160394
  idle_oscillation_freq_hz: 2596.0
  idle_oscillation_level_dbfs: -91.28972525799054
  transient_smear_ms: 0.0
  transient_overshoot_db: 2.4659084174875163
  transient_recovery_ms: 2.3125
```

## Extraction Progress Notes
- No significant ABI tracing bottlenecks noted in the source manifest.
- The module isolates correctly as a deterministic state machine block through direct signal analysis techniques, circumventing the need for raw binary extraction of individual CPU vectors. 

# Integration notes
- Implement the CVSD simulation entirely using a simple variable-step-size integrative delta block configured with these exact $\tau$ time constants.
- The signature sound fundamentally relies on `syllabic_charge_tau_ms` reacting slower than the `discharge` counterpart. Preserve this explicit 39ms/16ms lag inequality to accurately reconstruct the transient squashing heard in vintage samplers.
