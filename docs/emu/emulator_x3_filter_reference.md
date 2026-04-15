# Emulator X3 filter bank reference

Archived from the *Emulator X3 Reference Manual*, Chapter 6 "Voice
Processing — Amplifier, Filter & Auxiliary Envelopes." This is the
canonical filter bank TRENCH's baked presets come from. Do not
paraphrase — this is the ground truth.

## Filter bank overview

You can choose from 55 different filter types or choose "No Filter,"
which bypasses the filter section. Most of the filters have two
parameters: Frequency (or Morph) and Resonance (Q, Gain, Body Size).
These two parameters can be continuously varied during the note.

The frequency response curve of the filter is accurately shown in the
display as the initial filter settings are changed. Frequency is shown
on the horizontal axis, with amplitude on the vertical axis.

High values of Q or resonance amplify frequencies near the cutoff or
center frequency. In a swept EQ filter, gain controls the amount of
boost or cut. In a phaser or flanger, resonance determines the depth of
the effect. In the vocal filters, body size determines the apparent
size of the mouth cavity.

12th-order filters use more CPU and therefore decrease the maximum
voice count.

## Filter type codes

| Code | Description       |
|------|-------------------|
| BPF  | Band-pass         |
| DST  | Distortion        |
| EQ+  | EQ boost          |
| EQ-  | EQ cut            |
| FLG  | Flanger           |
| HPF  | High-pass         |
| LPF  | Low-pass          |
| PHA  | Phaser            |
| PROG | Programmable      |
| REZ  | Highly resonant   |
| SFX  | Special effect    |
| VOW  | Vowel / formant   |
| WAH  | Wah               |
| WOW  | Wah-Wah pedal     |

## 2nd / 4th / 6th-order filters (low-order bank)

| Name                        | Order | Type | Description                                                   |
|-----------------------------|-------|------|---------------------------------------------------------------|
| 2-pole Lowpass              | 02    | LPF  | Typical OB-type low-pass. 12 dB/octave slope.                 |
| 4-pole Lowpass              | 04    | LPF  | Classic analog synth low-pass. 24 dB/octave rolloff.          |
| 6-pole Lowpass              | 06    | LPF  | Steeper than 4-pole. 36 dB/octave rolloff.                    |
| 2-pole Highpass             | 02    | HPF  | 12 dB/octave slope.                                           |
| 4-pole Highpass             | 04    | HPF  | Classic 4-pole high-pass. Cutoff sweep cuts 4th-order HPF.    |
| 2-pole Bandpass             | 02    | BPF  | 6 dB/octave rolloff on each side of the passband.             |
| 4-pole Bandpass             | 04    | BPF  | 12 dB/octave rolloff on each side.                            |
| Contrary Band pass          | 06    | BPF  | Novel band-pass: peaks and dips midway in the range.          |
| SweptEQ 1 octave            | 06    | EQ+  | Parametric, ±24 dB, 1-octave bandwidth.                       |
| Swept EQ 2→1 octave         | 06    | EQ+  | ±24 dB; bandwidth 2 oct at low end, 1 oct at high end.        |
| Swept EQ 3→1 octave         | 06    | EQ+  | ±24 dB; bandwidth 3 oct at low end, 1 oct at high end.        |
| Phaser 1                    | 06    | PHA  | Phase-shifter comb filter. Q varies notch depth.              |
| Phaser 2                    | 06    | PHA  | Slightly different notch frequency spacing.                   |
| FlangerLite                 | 06    | FLG  | Three notches; Q increases flanging depth.                    |
| Vocal Aah-Ay-Eeh            | 06    | VOW  | Vowel sweep Aah → Ay → Ee. Q = mouth cavity size.             |
| Vocal Ooh-Aah               | 06    | VOW  | Vowel sweep Oo → Oh → Ah. Q = mouth cavity size.              |
| Dual EQ Morph               | 06    | PROG | 2-frame morphing filter, two EQ sections.                     |
| Dual EQ + Lowpass Morph     | 06    | PROG | Two EQ sections + LPF tied to Morph.                          |
| Dual EQ Morph + Expression  | 06    | PROG | Two EQ sections + independently-controllable LPF.             |
| Peak/Shelf Morph            | 06    | PROG | 2-frame morph, independent freq/shelf/peak per frame.         |
| Morph Designer              | 2-12  | PROG | Up to 6 programmable sections (LP/HP/EQ), 2 frames.           |

## 12th-order filters (high-order bank)

| Name          | Type | Description                                                          |
|---------------|------|----------------------------------------------------------------------|
| Ace of Bass   | EQ+  | Bass-boost to bass-cut morph.                                        |
| MegaSweepz    | LPF  | "Loud" LPF with a hard Q. Tweeters beware.                           |
| EarlyRizer    | LPF  | Classic analog sweeping with hot Q and low-end.                      |
| Millennium    | LPF  | Aggressive low-pass. Q gives a variety of spiky tonal peaks.         |
| MeatyGizmo    | REZ  | Filter inverts at mid-Q.                                             |
| KlubKlassik   | LPF  | Responsive LPF sweep with a wide spectrum of Q sounds.               |
| BassBox-303   | LPF  | Pumped-up lows with TB-like squelchy Q factor.                       |
| FuzziFace     | DST  | Nasty clipped distortion. Q = mid-frequency tone control.            |
| DeadRinger    | REZ  | Permanent "ringy" Q response. Many Q variations.                     |
| TB-OrNot-TB   | EQ+  | Great bassline "processor."                                          |
| Ooh-To-Eee    | VOW  | Oooh → Eeee formant morph.                                           |
| BolandBass    | EQ+  | Constant bass boost with mid-tone Q control.                         |
| MultiQVox     | VOW  | Multi-formant. Map Q to velocity.                                    |
| TalkingHedz   | VOW  | "Oui" morphing filter. Q adds peaks.                                 |
| ZoomPeaks     | REZ  | High-resonance nasal filter.                                         |
| DJAlkaline    | EQ+  | Band-accentuating filter. Q shifts "ring" frequency.                 |
| BassTracer    | EQ+  | Low Q boosts bass. Try sawtooth or square with Q = 115.              |
| RogueHertz    | EQ+  | Bass with mid-range boost and smooth Q. Sweep cutoff at Q = 127.     |
| RazorBlades   | EQ-  | Cuts a series of frequency bands. Q selects different bands.         |
| RadioCraze    | EQ-  | Band-limited "cheap radio" EQ.                                       |
| Eeh-To-Aah    | VOW  | "E" → "Ah" formant movement. Q accentuates peakiness.                |
| UbuOrator     | VOW  | Aah-Uuh vowel with no Q. Raise Q for throaty vocals.                 |
| DeepBouche    | VOW  | French vowels. "Ou-Est" vowel at low Q.                              |
| FreakShifta   | PHA  | Phasey movement. Try major 6 interval and maximum Q.                 |
| CruzPusher    | PHA  | Accentuates harmonics at high Q. Try with a sawtooth LFO.            |
| AngelzHairz   | FLG  | Smooth sweep flanger. Good with vox waves.                           |
| DreamWeava    | FLG  | Directional flanger. Poles shift down at low Q, up at high Q.        |
| AcidRavage    | REZ  | Great analog Q response. Wide tonal range. Try a sawtooth LFO.       |
| BassOMatic    | REZ  | Low boost for basslines. Q goes to distortion at the maximum level.  |
| LucifersQ     | REZ  | Violent mid-Q filter. Take care with Q values 40–90.                 |
| ToothComb     | REZ  | Highly resonant harmonic peaks shift in unison. Try mid Q.           |
| EarBender     | WAH  | Midway between wah and vowel. Strong mid-boost. Nasty at high Q.     |
| KlangKling    | SFX  | Ringing flange filter. Q "tunes" the ring frequency.                 |

## Morph Designer (PROG, 2–12 order)

The ultimate E-mu synthesizer filter. Up to six completely programmable
filter sections (lowpass, highpass, or EQ) and the ability to morph
between two completely different frame settings of the six sections.

### Architecture

```
    Hi Morph  [stage 1] [stage 2] [stage 3] [stage 4] [stage 5] [stage 6]
    Lo Morph  [stage 1] [stage 2] [stage 3] [stage 4] [stage 5] [stage 6]
                  LP        EQ        EQ        EQ        EQ        HP
                 (filter type per stage; user-selectable)
```

Any of the six filter sections can be designated as lowpass, highpass,
or EQ. Each section has a filter frequency and Q (or Gain in the case
of EQ) for both the Lo and Hi morph positions. When the Morph control
is changed, all the filter controls interpolate between the settings
programmed. Twenty-four knobs are reduced to two extremely powerful
controls: Morph and Gain (offset).

### EQ stage gain controls

- Small Gain knobs control Gain (−24 dB to +24 dB).
- The big Gain Wheel is Gain Offset (−24 dB to +24 dB).
- For any given section, the two Gain controls add (clipping at 0 or 100%).

### HP/LP stage gain controls

- Small Gain knobs control Q (0 to 100%).
- The big Q Wheel is Q Offset (−50% to +50%).
- For any given section, the two Q controls add (clipping at 0 or 100%).

### Control matrix

| Control             | Function                                                            |
|---------------------|---------------------------------------------------------------------|
| Morph (Frequency)   | Morphs between Lo Morph and Hi Morph frames. Modulation destination.|
| Gain/Q (Resonance)  | Adds to all section Gain/Q simultaneously. Modulation destination.  |
| Section             | Selects which filter section is currently being edited.             |
| Shape               | Off / EQ / Lowpass / Highpass for the current section.              |
| Lo Morph Frequency  | Center frequency at Morph minimum.                                  |
| Lo Morph Gain/Q     | Gain (EQ) or Q (LP/HP) at Morph minimum.                            |
| Hi Morph Frequency  | Center frequency at Morph maximum.                                  |
| Hi Morph Gain/Q     | Gain (EQ) or Q (LP/HP) at Morph maximum.                            |

## Other programmable filters

### Dual EQ Morph

A programmable 2-frame morphing filter with two EQ sections. As Morph
is increased, the filter interpolates from the Lo to Hi settings. The
Gain of each section remains constant during the Morph, but can be
scaled at note-on time by controlling the Initial Gain parameter. Can
be used to create custom vocal formant filters.

### Dual EQ Morph + Expression

Two EQ sections plus a separately controllable lowpass filter.
Expression controls the LPF frequency independently from Morph. Lets
velocity control expressive timbre via the LPF while Morph is modulated
for another effect. The gain parameter for each EQ section stays
constant during the Morph.

### 2EQ Morph + Lowpass Morph

Two EQ sections plus an LPF whose Fc is tied to Morph. At Morph = 255
the LPF is completely open. The LPF Q has an initial setting and can
also be modulated realtime via the Filter Resonance parameter. The gain
parameter for each EQ section stays constant during the morph.

### Peak/Shelf Morph

2-frame morph with independent control over frequency, shelving, and
peak for each frame. When Shelf is negative, low shelving response;
positive, high shelving; zero, peak response. Peak is a real-time
modulation destination via Filter Resonance.

## TRENCH implications

- A TRENCH pill is a baked Morph Designer preset: 6 sections × 2 frames
  × (section type, frequency, Q/Gain).
- The 30-integer heritage authoring grid at
  `cartridges/engine/_source/heritage_designer_sections.json` is the
  Morph Designer authoring format directly (6 sections × 5 fields).
- `TalkingHedz` is the factory 12-pole VOW filter baked into
  `trench-core/src/hedz_rom.rs` (soon `emu_params.rs`).
- `BassBox-303` is the factory 12-pole LPF listed here, NOT the
  audiorealism TB-303 emulator plugin.
- E-mu's Q/Gain is a live runtime parameter with a ±50% offset wheel,
  not two precomputed snapshots. TRENCH's current 4-corner shape
  precomputes Q0 and Q100 variants to avoid runtime biquad coefficient
  computation; see
  [`../architecture/zplane_truth.md`](../architecture/zplane_truth.md)
  for the 2-corner transition plan.
- 12th-order filters are legitimate (33 of the 55-filter bank), not
  artifacts. Stage count alone is not a reliable signal for extraction
  errors.
