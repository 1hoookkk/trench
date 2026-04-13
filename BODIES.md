# TRENCH bodies

The 4 launch bodies. Names are normative; everything else is the contract a
candidate must hit to ship under that name. **How to know it's wrong** lists
the falsifiable failure modes for each body — if any item triggers, the
candidate is rejected.

---

## 1. Speaker Knockerz

**Sonic identity:** highly pressurized sub-harmonic resonator that chokes
and fractures low-end weight into mid-range, tearing speaker-cone cry.

**Invariant:** the fundamental below 60 Hz remains anchored and phase-locked.
The sub never disappears, even as upper harmonics are brutalized.

**Motion logic:** slow creeping LFOs that occasionally snap via S&H into
"Cardboard Rip" territory.

**Differentiation:** treats bass as a physical air mass to squeeze, choke,
and tear — not as a frequency band to pass.

**Best source material:** saturated 808s, raw sine subs, plucky synth subs,
low male vocals.

**Morph path (notch territories):**
1. The Vault — clean sub, heavy weight, clamped highs.
2. Chest Resonance — sub blooms into a warm 150 Hz boost.
3. The Choke — steep notch at 200 Hz; fundamental separates from overtones.
4. Cardboard Rip — peaky unstable resonance flares at 400 Hz.
5. The Rattle — comb-like flutter in the 800 Hz range.
6. Cone Cry — narrow screaming peak at 1.2 kHz.
7. Total Fracture — mid-band phase-collapses; pure sub plus violent top.

**How to know it's wrong:**
- Sub below 60 Hz drops more than 6 dB at any (morph, Q) position.
- On a 50 Hz sine, fundamental level moves more than 3 dB across full morph.
- 200 Hz Choke notch is missing or shallower than 12 dB at its station.
- 1.2 kHz Cone Cry peak fails to exceed flanking bands by 10 dB.
- Raising Q produces broadband harshness instead of localized emphasis.
- Body sounds tonal rather than physical — no sense of pressure or tear.

---

## 2. Aluminum Siding

**Sonic identity:** brittle high-tension treble stressor that folds air and
sibilance into expensive crystalline damage.

**Invariant:** the midrange around 1 kHz remains scooped — a permanent
acoustic void. All kinetic energy lives in the extreme highs.

**Motion logic:** fast jittery envelope follow; snap into "Aluminium Tear"
on transients, release back to silk.

**Differentiation:** weaponizes the top end. Turns harshness and sibilance
into authored, musical texture.

**Best source material:** hi-hats, female vocals, overheads, granular pads,
foley.

**Morph path (notch territories):**
1. Dull Silver — muted phasey top with a strange airless void.
2. The Sheen — wide silky air boost around 10 kHz.
3. Glass Stress — narrow piercing peak at 7 kHz.
4. The Sibilant Fold — double-notch ~5 kHz and ~8 kHz folds S/T sounds.
5. Aluminium Tear — harsh metallic inharmonic ring at 12 kHz.
6. Dog Whistle — extreme resonance push at 18 kHz.
7. Shatter Point — upper boundary opens; burst of white noise/air.

**How to know it's wrong:**
- 1 kHz region rises above -12 dB relative to flanking bands at any (morph, Q).
- 12 kHz Aluminium Tear region sounds smooth/harmonic — no inharmonic ring.
- On hi-hats, transients soften instead of sharpening.
- On vocal "ess" content at high Q, 5–8 kHz double-notch fails to fold.
- Body sounds dark or warm at any morph position.
- 18 kHz peak is absent or below the 10 kHz Sheen level.

---

## 3. Small Talk Ah-Ee

**Sonic identity:** biomechanical vocal cavity that mutates from a deeply
relaxed open throat into a strangled digital scream.

**Invariant:** exactly two dominant formants are active at all times,
mimicking human vocal cord separation, even when surrounding frequencies
are mangled.

**Motion logic:** envelope-follow + LFO; opens throat on loud passages,
closes on tails.

**Differentiation:** simulates muscular throat tension. Wet, fleshy,
anatomical.

**Best source material:** lead vocals, brass stabs, supersaws, Reese basses.

**Morph path (notch territories):**
1. The Yawn — deep "Oo" vowel, massive cavity, hollow.
2. The Hum — closed-mouth "Mm," nasal buzz, trapped air.
3. Open Ah — wide intelligible "Ah," pushed forward.
4. The Bite — sharp "Ee" through narrow tube, 2.5 kHz presence.
5. The Gag — first formant collapses, strangled high squeal.
6. Digital Rasp — formants intermodulate like vocal fry / bitcrush.
7. The Shriek — wide-open distorted "Aaa," max throat tension.

**How to know it's wrong:**
- Sustained sawtooth shows fewer than 2 distinct resonant peaks at any
  (morph, Q) position.
- Formants drift outside 200 Hz – 4 kHz outside the Gag/Shriek extremes.
- Low-Q sweep sounds like a single-peak wah instead of two-formant motion.
- Open Ah on a 200 Hz saw doesn't read as "ah" by ear.
- Both formants move in lockstep instead of independently.
- Body sounds synthetic-EQ instead of throat-like.

---

## 4. Cul-De-Sac (the meta body)

**Sonic identity:** paradoxical structure that begins as a thick blunt
physical tube and physically fractures halfway through the morph into a
scattered multi-band comb matrix.

**Invariant:** a constant low-level resonant hum at the root note. The only
tether to reality while upper harmonics undergo structural collapse.

**Motion logic:** slow deliberate macro sweeps via mod wheel or unsynced
LFO. The user plays the State Boundary, teetering on the fracture.

**Differentiation:** two physical models (single thick resonance vs 8-peak
comb) connected by a chaotic phase-cancellation event. A filter that turns
into another filter.

**Best source material:** drum loops, full mix busses, neuro basses, dense
pads.

**Morph path (notch territories):**
1. Iron Pipe — rigid metallic tube, thick, blunt, localized.
2. The Rust — deep notches in lower mids thin the pipe.
3. The Bulge — unnatural swell at 3 kHz, structure bends.
4. State Boundary — total phase cancellation, sound vanishes.
5. The Fracture — audio reappears as 8 narrow comb peaks.
6. Glass Shards — comb peaks spread chaotic, non-harmonic.
7. Crystal Dust — peaks thin to slivers, granular shimmer.

**How to know it's wrong:**
- Root resonance moves with morph instead of staying fixed at the root.
- Iron Pipe (m=0) has more than one dominant peak.
- Glass Shards (m≈0.85) shows fewer than 6 distinct comb peaks.
- State Boundary (m≈0.5) does not produce an audible null on a sustained tone.
- The fracture is gradual instead of a sudden phase event.
- Comb peaks at Glass Shards are harmonically aligned with the root.

---

## Release gate (fail-closed)

Before shipping cartridge updates:

    cargo xtask release-gate

Rebuilds the emotion-stress corpus/splits, runs holdout promotion with
balanced sampling, and fails unless all 4 shipping bodies are
`holdout_promoted=True` and `holdout_gate_ok=True`.

Only the 4 names above are normative: `Speaker Knockerz`, `Aluminum Siding`,
`Small Talk Ah-Ee`, `Cul-De-Sac`. Measurable release logic lives in code
(`pyruntime.analysis.shipping_gate`), not prose.
