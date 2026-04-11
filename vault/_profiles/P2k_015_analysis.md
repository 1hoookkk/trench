## Flanger Lite (P2k_015, filterType=15)

### What It Is
This is one of the few names here that is not completely dishonest, but it is still not a classic jet sweep. It behaves like a bright comb chassis with one low throat bone underneath it. At low Q the upper stages lock into a hard 3-14 kHz teeth-and-air lattice; at high Q a 400-500 Hz floor appears while the rest of the structure stays metallic and carved. A producer would reach for this when a sound needs to feel scraped, hollow, and machine-tight, with pressure in the gaps rather than a big singing vowel.

### Stage Roles
- Stage 1: air lid at 11.6-15.0 kHz. Moves with Q. Stress.
- Stage 2: anchor at 427-1794 Hz, from AE/EH mouth down to vocal-tract floor. Moves with Q. Anchor.
- Stage 3: cavity/bridge at 3.1-6.6 kHz, crossing presence into sibilance. Both. Cavity.
- Stage 4: bite blade at 13.5-17.8 kHz, then slightly lower 15.4 kHz glass. Both. Stress.
- Stage 5: air bridge at 5.4-9.8 kHz, from upper bite to sibilant roof. Moves with Q. Bridge.
- Stage 6: comb hinge at 8.6-12.4 kHz, then 11.5 kHz against a far air zero. Both. Identity.

### Motion Law
Q is the real mechanism. At both morph extremes, raising Q drags Stage 2 down from 1.8 kHz to roughly 430-500 Hz while Stages 1 and 4 climb higher and Stages 3, 5, and 6 spread across the 6-12 kHz roof. That does not make the filter more vocal; it creates a low body under a sharpened metal grate. Morph at Q0 barely relocates anything. It mostly hardens the same anatomy by pushing all six radii toward `0.999`, which is why `M100_Q0` gets the biggest peak and the lowest centroid without changing the pole map much. Morph at high Q is subtler but more important for carve depth: the poles stay in roughly the same bands, yet the zeros retune hard enough to drive notch depth from about `-37 dB` to `-84 dB`. The most violent transition is `M100_Q0 -> M100_Q100`: the boosted comb body loses about `24 dB` of peak, dynamic range expands by about `43 dB`, and the whole surface turns from resonant grin into a much more surgical notch field. It feels physical because one throat bone stays low while the roof above it keeps tightening like metal teeth closing around a tube.

### What Would Break It
Stage 2 is the load-bearing piece. Without that 427-503 Hz floor at high Q, the filter stops feeling mounted to a cavity and becomes anonymous top-end combing.

### Character Tags
[metallic] [comb] [subtractive] [pressure] [morphing]

### Relatives
Closest to P2k_000 and P2k_002 because all three use a low anchor under a carved high-band shell, but `Flanger Lite` is the cleanest comb anatomy of the set: less speech, more teeth. Opposite of P2k_013, which spends its energy on coherent mouth formants; this one prefers a rigid grate with only a token throat underneath.
