## Peak/Shelf Morph (P2k_021, filterType=21)

### What It Is
This is an EQ-shaped “mouth box” that can act like a shelf/peak tilt or become a brutally hollow cavity with ribs of presence and near-silence between them. At `M0_Q0` it’s bright (peak ~8.4 dB) but carved deep (notches ~-98 dB); push morph and the centroid drops ~1.2 kHz while the voids reach ~-124.5 dB. Q is the twist: at low morph it vents into a smoother plate (range ~56 dB), but at high morph it snaps into a chest+throat stack under a 9.2 kHz spit. Use it when you want “boxed speech” that can turn surgical without becoming a static wah.

### Stage Roles
- Stage 1: bite/air lid at 3.6–9.2 kHz. Both. Stress.
- Stage 2: anchor at 123–580 Hz (chest ↔ F1). Both. Anchor.
- Stage 3: formant hinge at 216 Hz–2.14 kHz (throat floor ↔ F2 mask). Both. Identity.
- Stage 4: bridge from mask to sibilance at 1.95–9.32 kHz (trades places with Stage 1 at high Q). Both. Bridge.
- Stage 5: bite/presence plate at 1.57–3.92 kHz. Both. Stress.
- Stage 6: cavity/shelf knife at 0.85–2.84 kHz, backed by a near-unit zero that lives up at ~5–18 kHz. Both. Cavity.

### Motion Law
Morph mostly darkens and increases carve. At Q0 it drags Stage 2 from a real F1 (~580 Hz) into chest (~123 Hz) and pulls Stage 6 from ear-canal (~2.84 kHz) toward a mouth hinge (~850 Hz) while Stages 3/5 climb into a tighter 2–4 kHz mask, dropping centroid ~1169 Hz and deepening notches to ~-124.5 dB.

Q is a rewire switch. At `M0`, raising Q spreads poles into ≈200/890/1.57/2.35/4.61/9.32 kHz and vents the comb (dynamic range -50.5 dB). At `M100`, raising Q assembles ≈216/567/1.14/1.95/3.92 kHz plus a 9.2 kHz blade. The most violent move is `M0_Q100 -> M100_Q100`: Stage 4 falls 9.32 kHz → 1.95 kHz while Stage 1 jumps into sibilance, swapping where the “teeth” live.

### What Would Break It
Stage 6. Without that near-unit-circle zero that pins the shelf edge and makes the extreme voids possible, the filter reads like a busy midrange mask instead of a physical box with vents.

### Character Tags
[vocal] [surgical] [comb] [morphing] [pressure]

### Relatives
Closest to `P2k_018`–`P2k_020` (face assembly under Q) and a cousin to `P2k_015` for comb voids, but unique because Q at low morph *smooths* instead of tightening. Opposite of lid-like lowpass types (`P2k_001`–`P2k_003`).
