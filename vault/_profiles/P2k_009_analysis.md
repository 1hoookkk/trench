## Swept EQ 1 Octave (P2k_009, filterType=9)

### What It Is
This is not behaving like a polite swept bell. It feels like a throat cavity with a detachable metal roof: one state gives you a low spoken body, another strips that body away and leaves a hard air-band shell. A producer would reach for it when a sound needs to pivot between voiced pressure and brittle glare without sounding like a normal EQ sweep.

### Stage Roles
- Stage 1: throat pinch to air-band lid at 0-13.5 kHz, with a 357 Hz floor/zero memory. Both. Stress.
- Stage 2: mouth bridge at 1.1-1.9 kHz, then fractures into a 10.3 kHz roof at full morph/high Q. Both. Bridge.
- Stage 3: load-bearing formant from 357-376 Hz chest/throat up to 2.05 kHz mouth bite. Both. Identity.
- Stage 4: detached shell at 15.0-15.7 kHz, except one Q-raised collapse to DC with a 2.0 kHz cavity zero. Moves with Q. Cavity.
- Stage 5: nasal-to-mouth carrier at 753 Hz to 1.16 kHz, then jumps to a 13.3 kHz metallic crown. Both. Identity.
- Stage 6: anchor that swaps from 357-445 Hz throat floor to 1.97 kHz bite, then to 17.9 kHz air fracture. Both. Anchor.

### Motion Law
At `M0`, Q partly opens the filter upward: Stage 1 leaps from 578 Hz to 9.8 kHz and Stage 4 falls out of the audio band, but Stages 3, 5, and 6 stay in bodily territory around 445 Hz to 2.05 kHz. That gives you a mouth under a lifted roof. Morph at `Q0` does the opposite. Stage 3 drops from 1.29 kHz to 357 Hz while Stage 6 climbs from 357 Hz to 1.97 kHz, so the inner anatomy cross-swaps; meanwhile Stage 1 dies to DC and Stage 4 remains a 15 kHz shell. At `M100`, Q is the real fracture point: Stages 2, 5, and 6 rip from 1.1-2.0 kHz into 10.3-17.9 kHz while Stage 3 stays pinned near 376 Hz. The most violent transition is `M100_Q0 -> M100_Q100`, where a resonant throat-and-mouth stack loses body and becomes one low bone under a bright metal canopy. It feels physical because one part stays as the body while the rest peel upward like jaw and palate separating.

### What Would Break It
If Stage 3 weakened, the whole filter would lose its body reference. It is the one stage that stays load-bearing in low chest/throat territory even when the rest either cross-swap or detonate upward. Without it, the filter turns into top-end combing with no anatomy underneath.

### Character Tags
[metallic] [pressure] [fracture] [morphing] [biological]

### Relatives
Closest to P2k_008 and P2k_006 because all three split between voiced body and brighter shell, but P2k_009 is less orderly: it cross-swaps its low stages at low Q, then tears three stages into the air band at high morph/high Q. Opposite of P2k_007, which ratchets upward as a linked mouth; this one keeps threatening to separate the roof from the throat.
