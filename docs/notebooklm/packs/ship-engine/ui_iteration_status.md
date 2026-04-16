---
name: UI iteration status
description: Current state of plugin UI rewrite — response curve rendering, what works, what's broken
type: project
---

Plugin UI is being rebuilt as a frequency response curve display (480x320 landscape, wgpu).

**What works:**
- Standalone mode: `cargo run -p trench-plugin --bin trench-plugin --release`
- Background shader renders (dark void + grain + scanlines)
- Response curve mesh builds correctly (fill + trace + ghost layers)
- Curve shape is correct — peaks/cuts show in the right places
- Spectrum analyzer (FFT) wired into audio thread, lock-free AtomicU32 bridge
- Envelope follower always-on (env_active param removed)
- Scroll cycles bodies (guarded during drag)

**What's broken:**
- Fill shader is too dim — flat sections (0dB / full-cut) are invisible against the dark background
- Multiple shader iterations all produced near-identical output because fill color at p=0 is too close to background color
- Need to make a DRAMATIC visual change first to confirm shader changes are taking effect, then refine

**Key insight from user:**
- Stop making subtle incremental shader tweaks. Make something obnoxiously visible first, confirm it works, THEN refine.
- User wants to iterate one change at a time with visual validation each step.
- The response curve IS the right visualization. Just make it look good.
- Constellation/pole-zero overlay was a dead end ("trap remix youtube thumbnail")
- Game UI design vocabulary is the right framing (volumetric density, structural stress, vertex elevation)

**Why:** User has been frustrated by invisible curves and design tangents. The fill alpha-blends to near-zero on the dark background. Need to either: brighten the fill floor dramatically, or lighten the background in the display zone, or both.

**How to apply:** Next session should start by making the fill neon green at 80% opacity to prove the pipeline, then dial colors/opacity back to taste.

**Files:**
- `trench-plugin/src/ui.rs` — all rendering (shader, mesh, pipelines, window handler)
- `trench-plugin/src/constellation.rs` — pole-zero extraction (unused, `#![allow(dead_code)]`)
- `trench-plugin/src/spectrum.rs` — FFT analyzer (working, tested)
- `trench-plugin/src/response.rs` — frequency response computation (working)
- `trench-plugin/src/main.rs` — standalone entry point

**Release build takes ~2 min.** User gets impatient. Close old window BEFORE launching new build or exe lock blocks compile.
