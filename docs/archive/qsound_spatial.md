# Source
- `C:\Users\hooki\trench_re_vault\extract\qsound_spatial_v1\2026-03-05\`
- `C:\Users\hooki\trench_re_vault\extract\qsound_spatial_v1\2026-03-05_ingress_map_a\`
- `C:\Users\hooki\trench_re_vault\extract\qsound_spatial_v1\2026-03-06\`
- `C:\Users\hooki\trench_re_vault\extract\qsound_spatial_v1\2026-03-06_ingress_calc_a\`
- `C:\Users\hooki\trench_re_vault\datasets\qsound_spatial_v1\2026-03-05_pan_batch_fullgrid_a\`
- `C:\Users\hooki\trench_re_vault\models\qsound_spatial_v1\qsound_spatial_v1_recon_20260306_121748\`
- `C:\Users\hooki\trench_re_vault\datasets\qsound_spatial_v2_live_ingress\2026-03-06_smoke_b\`

# What this data is
This is tracking and model documentation for the reverse-engineered QSound spatial processing from Emulator X. It encapsulates offline `.wav` test grids, ABI failure traces for spatial endpoints, and the final deduced parametric clean-room math defining how az/el/distance translates into interaural time (ITD), level (ILD), and spectral shaping differences.

# Why it's here
QSound spatial placement is a first-class feature for the TRENCH architecture. Tokens placed into the Looperator grid won't just morph—they will have an explicit stereo/spatial pose. The CLEAN repo needs the derived ITD/ILD mathematical coefficients to seamlessly reimplement this 3D positioning natively.

# Sanitized content

## QSound Parameter Space
Based on ABI traces and parametric regression, the exposed QSound spatial topology accepts:
- **Pan**: -60 to +60
- **Distance**: Modeled via `log2(dist / 0.25)` curves.
- **Polar / Vector**: Azimuth (`az`), Elevation (`el`), and Distance.
- **Underlying Output Laws**: Synthesizes stereophonics through three correlated axes: ITD (microsecond ear delays), ILD (dB level skew), and 3-Band FIR Frequency Law filtering per ear. 

## V1 vs V2 Ingress Differences
- **V1 (Offline ABI)**: Relied on C# `.NET` interop targeting the `QMixer.dll` directly to batch render `.wav` files. However, while 1D pan setters worked, complex vector interfaces (`SetSourceVelocity`, `SetPolarPosition`, `SetSourcePosition`) threw fatal access violations (`0xC0000005`) internally within the host `QMixer.dll`.
- **V2 (Live Ingress)**: Pivoted entirely toward active live-plugin stimulus capture (`smoke` and `regression` tests) bypassing the broken standalone `.dll` vector host initializations entirely.

## Pan-Batch Measurement Grid
`pan_batch_fullgrid_a` measures simple linear array placement across exactly 13 vectors running cleanly from Left -60 to Right +60 in steps of 10.

| Value | Filename Capture | Role / Side |
|---|---|---|
| -60 | `pan_-060.wav` | Hard Left |
| -50 | `pan_-050.wav` | |
| -40 | `pan_-040.wav` | |
| -30 | `pan_-030.wav` | |
| -20 | `pan_-020.wav` | |
| -10 | `pan_-010.wav` | |
| 0 | `pan_+000.wav` | Center |
| +10 | `pan_+010.wav` | |
| +20 | `pan_+020.wav` | |
| +30 | `pan_+030.wav` | |
| +40 | `pan_+040.wav` | |
| +50 | `pan_+050.wav` | |
| +60 | `pan_+060.wav` | Hard Right |

## Canonical Recon Model
The definitive, latest clean-room solver run is **`qsound_spatial_v1_recon_20260306_121748`**.

### Target Coefficient Tables (Parametric Model)

```json
{
  "itd_law": {
    "feature_order": [
      "sin(az)", "sin(2*az)", "sin(3*az)", "sin(4*az)", "sin(5*az)", "sin(az)*abs(el)/30"
    ],
    "coef": [
      3578.7646232504208,
      -99.11300089026597,
      -960.7096166779854,
       631.0360877559084,
      -229.78233226297414,
      -173.11867492761178
    ]
  },
  "ild_law": {
    "feature_order": [
      "sin(az)", "sin(2*az)", "sin(3*az)", "sin(4*az)", "sin(5*az)", "sin(az)*abs(el)/30"
    ],
    "coef": [
      6.819731516349544,
      -2.5008130981129657,
      0.8210199212077876,
      -0.198121089137829,
      0.021401263894609834,
      8.819077411548193e-05
    ]
  },
  "band_law": {
    "feature_order": [
      "1", "log2(dist/0.25)", "log2(dist/0.25)^2", 
      "sin(az)", "cos(az)", "sin(2*az)", "cos(2*az)", 
      "sin(3*az)", "cos(3*az)", "el/30", "sin(az)*el/30", "cos(az)*el/30"
    ],
    "channels": {
      "l": {
        "low": [ -72.211726, 0.386310, -1.609965, -2.942365, 4.905630, 0.894425, -1.686776, -0.082666, 0.833466, -5.542e-17, 2.501e-16, 1.788e-16 ],
        "mid": [ -74.892487, 0.460364, -1.645200, -2.794224, 9.547364, 1.100615, -2.989885, -0.168168, 1.485241, 1.372e-17, 2.241e-17, 3.424e-16 ],
        "high": [ -73.911507, 0.438043, -1.634712, -2.835356, 8.046734, 1.041892, -2.546882, -0.132837, 1.302472, -3.615e-17, 3.924e-17, -1.875e-16 ]
      },
      "r": {
        "low": [ -72.211726, 0.386310, -1.609965, 2.942365, 4.905630, -0.894425, -1.686776, 0.082666, 0.833466, 6.880e-17, 9.784e-17, -3.152e-16 ],
        "mid": [ -74.892487, 0.460364, -1.645200, 2.794224, 9.547364, -1.100615, -2.989885, 0.168168, 1.485241, -3.922e-16, 1.316e-16, 9.925e-17 ],
        "high": [ -73.911507, 0.438043, -1.634712, 2.835356, 8.046734, -1.041892, -2.546882, 0.132837, 1.302472, -1.004e-16, 1.863e-16, -2.686e-16 ]
      }
    }
  }
}
```
*(Notes: ITD and ILD are scaled mathematically by a right-ear-lead/right-ear-louder convention. High-frequency target gains are computed strictly from the modeled 3-Band crossfeed coefficient matrix.)*

# Integration notes
- Looperator tokens should carry an optional `[az, el, dist]` spatial pose array attribute.
- The runtime token compiler will feed this `[az, el, dist]` pose into the canonical V1 model coefficients.
- Apply the modeled `band_law` as a 3-way multi-band FIR block post-filter, adjusting explicit L/R channel gain offsets per the `ild_law`, and introducing a relative sample delay across channels per the `itd_law`.

# UNKNOWNs
- None strictly blocking. The exact DSP topology for the 3-Band FIR wasn't embedded here, but the data cleanly dictates how pan maps into spatial properties at a token level.
