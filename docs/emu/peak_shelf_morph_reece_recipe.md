# Peak/Shelf Morph — Reece Stab Recipe

Executable recipe for the DnB reece stab configuration described in the
DillusionMan tutorial. Companion to
`docs/emu/dillusionman_peak_shelf_morph.md` (which is the full
vocabulary reference). Original text:
`trenchwork_clean/how to configure a peakshelf morph.txt`.

## Step 1 — Left Morph Frame (starting position)

Establishes the foundational dark tone. Keeps the opening of the filter
low while leaving headroom for saturation.

- `FREQ`  = `246 Hz`
- `SHELF` = `-50` (slightly above full `-64` low-pass to prevent pops
  and nasty distortion when the filter is heavily saturated)
- `PEAK`  = `-24 dB`

## Step 2 — Right Morph Frame (target position)

Dictates the frequencies that emerge as the filter opens. Creates the
aggressive mid-range character.

- `FREQ`  = `4488 Hz`
- `SHELF` = `+30`
- `PEAK`  = `+1.5 dB`

## Step 3 — Program the morph sweep

Map a controller to transition between the two frames via the sampler's
Cords matrix.

- `MidiA → FilFreq @ +100%`

With the Peak/Shelf filter selected, `FilFreq` controls the `MORPH`
value between frames. Rotating the CC controller to a higher value
produces a sudden belching sound — slightly high-passed and
mid-shelved — modulating out of the foundational low-passed tone.

## Step 4 — Saturate the filter

Drive the filter's resonance to give the reece its tearing,
harmonic-rich quality during the transition.

- `MidiB → FilRes @ +100%`

Pushing this to a higher CC value saturates the filter and adds extra
harmonics to the big sweeps.

---

*This is the reece stab worked example. For the full SHELF/FREQ/PEAK
vocabulary and the why behind each setting, see
`docs/emu/dillusionman_peak_shelf_morph.md`.*
