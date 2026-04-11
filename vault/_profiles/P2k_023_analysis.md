## Ace of Bass (P2k_023, filterType=23)

### What It Is
Despite the name, this isn’t a sub-bass “bigger lows” filter — the lowest pole it ever plants is ~0.75 kHz. It’s a bassline *articulation rig*: a tight throat/mouth band (0.75–2.3 kHz) bolted under a serrated roof (3–15 kHz) that can go from glossy mask to near-silence comb voids. When you want a bass to *speak* and cut (acid syllables, metallic bite, hollow pressure), this is the one.

### Stage Roles
- Stage 1: air roof at 9.77–15.65 kHz (air band lid that drops hard off the origin). Both. Stress.
- Stage 2: anchor throat bone at 0.75–0.87 kHz (upper F1 / “jaw” ring). Moves with Q. Anchor.
- Stage 3: formant carrier at 1.24–2.28 kHz (F2 mouth / tongue position). Both. Identity.
- Stage 4: air brace at 10.51–13.05 kHz (air band spar with a sibilance-facing zero at high Q). Both. Stress.
- Stage 5: bite→air slider at 3.24–10.32 kHz (presence peak ↔ air). Both. Stress.
- Stage 6: cavity tooth at 6.73–12.64 kHz, with a near-unit-circle zero that jumps ~8.98 ↔ 17.96 kHz. Moves with Q. Cavity.

### Motion Law
Think “throat + tongue + metal grille.” Stages 2–3 are the body (0.75–2.3 kHz); Stages 1/4/6 are the roof; Stage 5 is the hinge that decides whether the roof lives in presence or pure air.

At low morph, raising Q doesn’t just tighten — it *pulls the roof down* (Stage 1: 15.6 → 9.8 kHz; Stage 5: 10.3 → 4.4 kHz; Stage 6: 12.6 → 8.0 kHz) while the tongue lifts (Stage 3: 1.86 → 2.28 kHz). That counter-motion is why it reads as a moving mouth inside a metallic helmet, not a static EQ.

Morph is the fracture lever at Q0. Sliding toward `M100_Q0` explodes the carve (peak ~12.7 → 31.5 dB; notches ~-42.2 → -100.7 dB; range +77 dB): the presence hinge (Stage 5) drops into the mid-roof and Stage 1 collapses from ultrasonic air into ~9.8 kHz, as if the lid slams shut and the grille starts ringing. The most violent “snap back” is raising Q at high morph (`M100_Q0 -> M100_Q100`), where the whole grill vents and the dynamic range collapses ~49 dB.

### What Would Break It
Stage 6. Without that near-unit-circle tooth (and its big Q-driven relocation), the filter loses the deep voids and stops feeling like pressurized metal — it becomes a generic midrange formant with some brightness.

### Character Tags
[metallic] [comb] [pressure] [fracture] [morphing]

### Relatives
Closest in *texture* to `P2k_015` (high-band teeth over a single body bone), and closest in *surface behavior* to `P2k_022` (a comb-mask that reconfigures instead of “opening/closing”). Opposite of the lowpass-family bodies (`P2k_001`–`P2k_003`), which read like a roof sliding down over a real chest floor rather than a mid-throat helmet.

