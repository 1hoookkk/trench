# E-mu Z-plane Morph Designer — HTML Prototype

Pure HTML/JS prototype of the E-mu Z-plane Morph Designer for the TRENCH project.
Stdlib only: ES modules, Canvas2D, AudioContext, hand-written radix-2 FFT.

## Run

```
python -m http.server 8765 --bind 127.0.0.1
```

Then open: `http://127.0.0.1:8765/demo/emu_zplane/`

## What it does

- Loads all 83 heritage designer pills from `data/heritage_designer_sections.json`
- Compiles pills at runtime through the E-mu type 1/2/3 firmware recipes (ported from `tools/bake_hedz_const.py`)
- Runs compiled kernels through a 12-stage DF2T cascade (ported from `trench-core/src/cascade.rs`)
- Draws frequency response (log-freq, dB) and impulse waveform via Canvas2D
- Plays audio through the Web Audio API

## Controls

| Control | Description |
|---------|-------------|
| Pill | Select heritage filter template (83 available) |
| Morph | 0 = low endpoint, 1 = high endpoint |
| Q | Frequency shaping: 0.5 = neutral, 0 = -50% freq, 1 = +50% freq |
| Impulse/Sine/Noise/Click | Stimulus type for frequency analysis and playback |
| Play / Stop | Trigger/stop looped audio playback |
| Parity Check | Runs 4-corner cross-language parity test against golden data |

## Acceptance criteria

- Parity check: all 4 corners green (maxAbsErr < 1e-4)
- Sliders responsive (debounce 30 ms)
- No console errors on load

## Files

```
demo/emu_zplane/
  index.html          — entry point
  style.css           — Bassbox-303 aesthetic
  app.js              — main wiring
  cascade.js          — 12-stage DF2T cascade
  emu_compile.js      — E-mu firmware type 1/2/3 recipes
  interpolate.js      — morph/Q interpolation + compile
  fft.js              — radix-2 FFT (N=1024)
  render_audio.js     — stimulus generators + offline render
  plot.js             — Canvas2D frequency + impulse plots
  parity.js           — cross-language parity check
  pills.js            — load + validate heritage_designer_sections.json
  data/
    heritage_designer_sections.json  — 83 E-mu templates
    hedz_golden.json                 — 4-corner golden impulse responses
```
