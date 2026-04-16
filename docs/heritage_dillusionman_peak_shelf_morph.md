# DillusionMan Peak/Shelf Morph Guide (October 2005)

Source: DillusionMan, 22 October 2005. Community guide for E-mu sampler Z-Plane Peak/Shelf Morph filter.

## Key Parameters

### MORPH (FilFreq CC sweep)
- Starting position of cutoff between two filter frames
- Left frame = 0, Right frame = 255
- FilFreq CC sweeps between frames (0 = left frame, 127 = right frame)
- Setting one patch to 0 and duplicate to 255 with same CC = inverted filter layer

### SHELF (-64 to +63)
- Defines filter tone character for both morph frames
- -64 = 100% low pass
- 0 = mid shelf
- +63 = high pass
- Values between blend the types:
  - -32 = 50% low pass + 50% mid shelf
  - -45 = 75% low pass + 25% mid shelf
  - +32 = 50% mid shelf + 50% high pass

### FREQ
- Context-dependent based on SHELF setting:
  - Low pass (SHELF negative): defines rolloff frequency
  - Mid shelf (SHELF ~0): defines the band center
  - High pass (SHELF positive): defines boundary frequency
- "It acts as a controller for the filter type"

### PEAK
- Master filter volume
- Adjusts relative volumes between frames for the sweep
- Per-frame Peak controls the balance during morphing

## Example: DnB Reece Stab

| Parameter | Frame 1 (low morph) | Frame 2 (high morph) |
|-----------|---------------------|----------------------|
| Freq      | 246 Hz              | 4488 Hz              |
| Shelf     | -50                 | +30                  |
| Peak      | -24 dB              | +1.5 dB              |

- Frame 1 Shelf = -50: keeps a low tone at the filter open, set above -64 to avoid pops and distortion during big CC sweeps
- Frame 2 Shelf = +30 with Peak +1.5 dB: CC rotation produces "sudden belching sound, slightly highpassed and slightly mid shelved, modulating out of the existing low passed tone"

## MIDI Setup

```
Cord: MidiA → FilFreq +100%
Cord: MidiB → FilRes  +100%
```

- Check CC assignments in Master > Midi > Cntrls2
- FilFreq sweeps the MORPH between frames
- FilRes adds harmonics (saturation)

## Context

- Most prominent in pre-2000 DnB: Optical, Grooverider, Dillinja
- The "belching" bass sound = saturated Peak/Shelf Morph sweep
- Z-Plane morphing = smooth transition between two authored filter states
- This is the exact authoring paradigm TRENCH's morph axis inherits
