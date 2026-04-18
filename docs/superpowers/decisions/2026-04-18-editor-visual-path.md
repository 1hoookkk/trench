# Editor Visual Path — Decision Primer

**Question to decide:** Which rendering approach for the TRENCH v1 plugin editor
maximizes unit sales?

- **A.** Raster faceplate (use existing PNG art as fixed background, overlay vizia
  widgets). Half a day. Flat 2D.
- **B.** Custom-drawn vizia Canvas (hand-painted elements: rollers, LCD,
  faceplate panels). 1-2 weeks. 2D animated.
- **C.** Raw wgpu editor (shader-based, depth/lighting/glow/particles possible).
  2-3 weeks. Full 3D/motion.

Pick one. Justify by commercial logic, not aesthetic preference. Ship path is
already Rust (nih-plug + `trench-live`) — the editor framework is separate from
that call.

---

## What TRENCH is

A Z-plane filter plugin. Reverse-engineered E-mu Emulator X filter topology
(12-stage DF2T biquad cascade, 4-corner morph interpolation) with heritage
E-mu AGC table verified against `EmulatorX.dll`. Targets an "intricate, hostile
filter instrument" aesthetic — not a multi-FX rack.

Shipping as Windows standalone + VST3 + CLAP. Plans for macOS later.
v1 ship target: **2026-07-15**.

**Signal chain v1:**
```
input → pre-cascade Mackie-1202-modeled softclip (per-body authored drive)
      → 12-stage DF2T cascade (MORPH + Q, REV/ENV modulation)
      → QSound-modeled spatial (user SPACE knob 0..1)
      → heritage E-mu AGC post-cascade
      → output
```

**User-facing controls (final v1 face):**
- BODY selector (4 bodies: Aluminum Siding, Speaker Knockerz, Small Talk Ah-Ee,
  Cul-De-Sac)
- MORPH, Q
- SPACE, TRIG
- REV, ENV (buttons)

Drive is baked per body. No DRIVE knob. No distortion block. No mod matrix.

## Who buys TRENCH

**Primary audience:** Underground trap / drum-and-bass producers chasing
"Reese stabs that feel physical." The producer-memory persona:
- Spends $30–$200 per plugin on boutique titles
- Buys based on YouTube demo videos (visual hype + audio demo)
- Shares plugin finds on Twitter/Discord — screenshots are the hook
- Already owns Serum/Vital (wavetable), Portal (granular), Phase Plant
  (modular), Output Arcade (loop), Pigments (hybrid), Diva (analog emu).
  TRENCH must earn its slot against these, not slot in behind.

**Secondary audience:** Sound designers and heritage-E-mu synth enthusiasts.
Smaller group, higher willingness-to-pay.

## Where we are right now

- **Rust ship path proven.** `target/release/trench-live.exe`, `trench-live.vst3`,
  `trench-live.clap` all build and launch. Editor window opens, event loop runs,
  MORPH/Q sliders wire via `ParamSlider` (nih-plug's canonical pattern).
- **DSP:** cascade + AGC + function_generator shipped in `trench-core/`.
  75 tests passing. Mackie drive and QSound spatial blocked on separate
  research issues (model source, measurement data quality).
- **Editor currently shows:** vizia defaults. Grey rectangles as body buttons.
  Mint rectangle labeled "LCD". Default `ParamSlider` widgets. Red and grey
  button placeholders. User feedback verbatim: **"looks generic as shit."**

## The three paths in detail

### Path A — Raster faceplate background + overlay vizia widgets

`assets/faceplate_blue.png` already exists (authored iteratively by a prior
Codex session). The current JUCE plugin uses it as a fixed background, with
roller components drawn over specific bounds.

Port that same pattern to vizia: load the PNG as a background image, overlay
`ParamSlider` (or custom-styled) widgets at matching coordinates.

- **Time:** ~half a day.
- **Identity match:** 100% — same art the JUCE plugin shipped.
- **Motion:** static PNG background. Rollers move as sliders. LCD is a cropped
  region of the image (or a separate widget overlay).
- **Demo video impact:** average. Screenshots work. Video feels 2005-era.
- **Differentiator strength:** low. Lots of plugins use this pattern.
- **Upgrade path:** drop-in replace with B or C later.

### Path B — Custom-drawn vizia Canvas

Build every visible element by hand using vizia's Canvas API. Rollers drawn
as rotating wheels with tick marks. LCD draws a live spectrum trace or
scanning grid. Faceplate panel drawn procedurally with gradient/shadow
treatment. Stencil font for labels. Visible screws, PCB hints.

nih-plug's own **Spectral Compressor** uses vizia Canvas for its complex
spectrum display — so this is a proven pattern in the ecosystem.

- **Time:** 1-2 weeks (mostly art direction + iteration, not code).
- **Identity match:** high if the art direction lands.
- **Motion:** animated elements, smooth rolling wheels, responsive LCD.
- **Demo video impact:** strong. In-motion GIFs land. Feels "alive."
- **Differentiator strength:** medium-high. Few boutique plugins do this
  level of hand-drawn 2D; most fall back to raster assets + standard widgets.
- **Upgrade path:** Canvas work informs wgpu work if Path C happens later.

### Path C — Raw wgpu editor (shader-based)

Full control rendering via wgpu + baseview. Custom shaders for the LCD glow,
plastic-case specular highlights, LED emission, particle burst on note hits,
depth/lighting for an actual 3D faceplate. Use the `wgpu-rendering` skill.

- **Time:** 2-3 weeks for the shader pipeline + editor; additional iteration
  for tuning.
- **Identity match:** max — can do things raster and Canvas can't.
- **Motion:** full procedural animation, shader-driven, 60fps buttery.
- **Demo video impact:** highest. Serum 2 / Phase Plant / Arcade tier of
  visual polish.
- **Differentiator strength:** very high. Very few plugins ship with
  purpose-built wgpu editors in 2026.
- **Upgrade path:** the end state. Nothing beyond.

## Commercial benchmarks (2026 market context)

Plugins that sell on visuals alone:
- **Serum** — default-grade UI but branded, ubiquitous. Sells via reputation.
- **Phase Plant** (Kilohearts) — bold flat design, sold by motion demos.
- **Portal** (Output) — pure visual hook, grew a fanbase on Instagram clips.
- **Arcade** (Output) — subscription model, heavy video marketing.
- **Pigments** (Arturia) — slick 3D modules, strong marketing budget behind.
- **Diva** (u-he) — classic analog skeuomorphism, sold via audio quality more
  than visuals.
- **Portal, Pattern, Thermal** (Output family) — very strong distinctive
  identities that define the product.

Plugins that sell despite default visuals:
- **Soundtoys Decapitator** — 2000s-era skeuomorphism, kept selling on sound
  quality and reputation.
- **FabFilter** line — clean/technical, sells on UX clarity not visual bling.

Plugins that failed despite good sound:
- Countless "boutique" plugins with generic default widgets — don't get past
  the "first impression" wall on YouTube thumbnails.

## Marketing vehicle matters

- Trap YouTube demo channels (Andrew Huang, Venus Theory, Cameron Clark, etc.)
  show **plugin in motion** while audio plays. Static faceplates get less
  engagement.
- Discord/Twitter audio-community shares rely on **short video clips, 15-30
  seconds**. Motion identity sells.
- Product landing pages feature a **hero GIF** — top 10% of plugins have
  motion that sells the "holy shit" moment visually.

## Skills available to the implementer

- `juce-lookandfeel-design` (not applicable — we're on Rust, not JUCE)
- `wgpu-rendering` (direct support for Path C)
- `generate-image` (AI image generation — can generate faceplate art iterations)
- `scientific-schematics` (for any technical/diagrammatic elements)
- Standard Rust + vizia + nih-plug tooling

## Constraints

- **Ship date:** 2026-07-15. ~3 months out.
- **Solo operator:** the user writes zero code. Subagents execute.
- **Single source of truth doctrine:** one codebase, one build, no parallel
  implementations.
- **Doctrine on taste:** no scorecards on user's creative calls; the user
  picks art by eye. Implementation decisions (framework, rendering approach)
  are engineering calls that must optimize for the instrument succeeding.

## Your call

Pick A, B, or C. Answer with:
1. **Pick.** Single letter.
2. **Why it sells the most units.** 3-5 bullet points, commercial reasoning
   only. No aesthetic preference language.
3. **Risk it fails.** The one failure mode that would kill this pick.
4. **What the first week of work looks like.** Concrete deliverables per day.
5. **Deferred work.** What gets pushed past v1 because of this pick.
