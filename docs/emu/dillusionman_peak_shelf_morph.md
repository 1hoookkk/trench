# EMU Peak/Shelf Morph — DillusionMan Tutorial

Heritage authoring vocabulary for the E-mu Peak/Shelf Morph filter, as
explained by DillusionMan (22 October 2005). This is the canonical
Z-plane authoring language TRENCH inherits. Source PDF:
`EMU Dillusion_PeakShelfMorph_Tutorial_WEB.pdf` at repo root.

## Screen reference

Two-frame filter, Peak/Shelf Morph. Example settings shown in the PDF:

| Control | Low Morph frame | High Morph frame |
|---|---|---|
| Freq  | 246 Hz | 4488 Hz |
| Shelf | -50    | +30    |
| Peak  | -24 dB | +1.5 dB |

Initial offset: Peak `-24 dB`, Morph `0`.

## The five controls

### 1. MORPH
Starting position of cutoff (FilFreq) value 0 within the two filter
frames.

- Left frame  = 0
- Right frame = 255

Setting either extreme provides the most drastic filter. Duplicate a
patch, set one to 0 and the other to 255 with the same CC modulations,
and you have inverted the filter and layered it with the original.

### 2. PEAK (master)
Master filter volume. Altering this modifies the peak filter's volume
as a whole.

### 3. FREQ
No fixed definition. Its effect depends entirely on the `SHELF` value:

- For a low-pass setting, `FREQ` defines the low-pass rolloff.
- For a mid-shelf setting, `FREQ` defines the band.
- For a high-pass setting, `FREQ` defines the high-pass boundary.

`FREQ` acts as a controller for the currently selected filter type.

### 4. PEAK (per-frame)
Adjust this value within each of the two frames to set the relative
volumes between them for the filter sweep.

### 5. SHELF
Defines the filter tone for each morphing frame. The figure represents
the range of the filter:

- `-64` = low pass
- `  0` = mid shelf
- `+63` = high pass

Values between these boundaries blend the aspects of the filter types.

Approximate blend examples:

| Value | Blend |
|---|---|
| `-64` | 100% low pass |
| `-45` | 75% low pass, 25% mid shelf |
| `-32` | 50% low pass, 50% mid shelf (emphasis on LP) |
| `-20` | mostly mid shelf with some low pass |
| `  0` | 100% mid shelf |
| `+32` | 50% mid shelf, 50% high pass |
| `+45` | 25% mid shelf, 75% high pass |
| `+63` | 100% high pass |

## Uses of FREQ and SHELF

- **Low-pass a sub:** `FREQ = 63 Hz`, `SHELF = -63`.
- **High-pass a lead / percussion:** `FREQ = 6000 Hz+`, `SHELF = +63`.
- **Cut mid from a lead** that's taking weight from the snare:
  `FREQ = 1000–4000 Hz`, `SHELF = 0`.

## Example — pre-2000 DnB reece stab

From the screen reference above:

- 1st frame `SHELF = -50` keeps a low tone at the open of the filter.
  It is deliberately not at `-64` so it avoids pops and nasty
  distortions when the filter is saturated by big CC sweeps.
- 2nd frame `SHELF = +30` with `PEAK = +1.5 dB` means a CC controller
  rotated to a higher value produces a sudden belching sound, slightly
  high-passed and slightly mid-shelved, modulating out of the existing
  low-passed tone.

Artists exploiting this sound in that era: Optical, Grooverider,
Dillinja.

## Cord setup

Best modulation uses the `MidiA`/`MidiL` controllers. Check what CC
values these are assigned to in `Master > Midi > Cntrls2` and map them
on your MIDI controller.

Cord layout:

```
MidiA   FilFreq   +100%
MidiB   FilRes    +100%
```

## Usage

With Peak/Shelf selected, `FilFreq` sweeps the `MORPH` value between
the two frames. With `Morph = 0`:

- `FilFreq = 0`   → left frame
- `FilFreq = 127` → right frame

Moving between these CC values creates the Z-plane effect of morphing
the two filter frames into a smooth transition.

## Saturating the filter

Set `FilRes` to a higher CC value to add harmonics. `FilRes = 127`
gives the most drastic filter — it acts like a Q on the overall sweep.

---

*Transcribed and reorganized from DillusionMan's infographic tutorial
(22 October 2005). Treat this as heritage authoring vocabulary, not
shipping TRENCH math. See `docs/emu/zplane_explained.md` and
`docs/architecture/zplane_truth.md` for how TRENCH relates to the
E-mu native model.*
