## 4 Pole Lowpass (P2k_002, filterType=2)

### What It Is
This name is also misleading. It is not behaving like a straightforward lowpass; it behaves like a combed throat that can flip between two anatomies: a bright ultrasonic tooth grid with only a low chest murmur underneath, or a stacked 260 Hz to 3 kHz vocal tract with deep carved nulls above it. A producer would reach for this when they want a sound to feel spoken through bone and metal, especially when sweeping from hollow mutter into a harder, more intelligible mouth.

### Stage Roles
- Stage 1: air at 14.9-15.5 kHz. Moves with Q via its zero; morph barely matters. Stress.
- Stage 2: anchor at 270 Hz or air at 15.3 kHz. Moves with Q. Identity.
- Stage 3: bite at 2.5-2.9 kHz or sibilance at 9.6-9.7 kHz. Moves with Q. Stress.
- Stage 4: bridge at 2.1-2.2 kHz or air at 13.2-13.4 kHz. Moves with Q. Identity.
- Stage 5: formant at 393-1236 Hz, from chest/OO floor into AH territory. Both. Anchor.
- Stage 6: anchor at 264-313 Hz or air at 10.9-11.8 kHz. Moves with Q. Load-bearing anchor.

### Motion Law
Q is the hinge. At Q0, Stages 1, 3, 4, and 6 sit high in the sibilance/air band while only Stages 2 and 5 touch the body, so the filter reads like a brittle metal shell with a faint low throat. Raising Q does not simply tighten resonance; it inverts the anatomy. Stage 2 drops from 15.3 kHz to 270-313 Hz, Stage 6 falls from 11-12 kHz to 264-313 Hz, and Stages 3-4 collapse together into a 2.0-3.0 kHz mouth/presence cluster. Morph then decides whether that body whispers or pushes. At Q0, morph drags the centroid down by 1.1 kHz and creates the only real boost state, because Stage 5 rises from 393 Hz to 980 Hz while the low anchor stays put. At Q100, morph keeps the same overall anatomy but shifts the floor: Stage 1/Stage 5/Stage 6 swap the sub-throat weighting enough to lift the centroid back upward. The most violent transition is M100_Q0 to M100_Q100, where the 9.5 dB boosted shell loses its grin and reassembles as a fully voiced cavity. It feels physical because the filter literally trades teeth for throat.

### What Would Break It
If Stage 6 weakened, the filter would lose the body-flip that makes the whole structure read as a mouth instead of a bright comb. Stage 6 is the load-bearing piece because Q uses it as the second low anchor under the 2-3 kHz mouth.

### Character Tags
[vocal] [metallic] [comb] [morphing] [pressure]

### Relatives
It is closest to P2k_001 because both use a low anchor under a metallic upper lid, but P2k_002 is more binary: Q swaps entire stages between air and throat instead of merely pulling them downward. P2k_000 is the opposite emphasis so far: more shredded air fracture, less coherent chest-and-mouth rebuild.
