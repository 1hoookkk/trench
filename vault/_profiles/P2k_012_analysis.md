## Phaser 1 (P2k_012, filterType=12)

### What It Is
This is not a classic phaser sweep. It behaves like a split mouth-and-metal chassis that can either keep its body low and let the roof hiss above it, or snap into a harder mask where one chest note, one upper-mouth note, and a high blade fight for control. A producer would reach for it when a sound needs to feel like air moving through bone and sheet metal at once: less swirl, more talking pressure with comb scars.

### Stage Roles
- Stage 1: air blade at 6.3-10.3 kHz. Both. Stress.
- Stage 2: bridge that flips from 274 Hz chest to 4.8 kHz bite. Both. Bridge.
- Stage 3: presence/formant spine at 2.9-3.5 kHz. Stays put. Identity.
- Stage 4: bite roof at 4.5-6.2 kHz. Both. Stress.
- Stage 5: identity hinge from 1.8 kHz upper mouth to 8.2-9.1 kHz air lid. Moves with morph. Identity.
- Stage 6: anchor from 225-412 Hz chest floor to 933-1103 Hz throat. Both. Anchor.

### Motion Law
The core relationship is a cross-swap, not a sweep. Stage 3 stays nailed near 3 kHz, acting like a fixed faceplate. Around that plate, morph trades jobs between Stage 5 and the upper cluster: at low morph, Stage 5 lives in the 8-9 kHz air band while Stage 6 holds the only real body at 225-1103 Hz; at high morph, Stage 5 drops into the 1.8 kHz mouth while Stages 1 and 4 harden into a brighter 6-10 kHz shell. Q then decides whether Stage 2 joins the body or the blade. At `M0`, raising Q drags Stage 2 down from 1.8 kHz to 274 Hz and lifts Stage 6 from chest to throat, so the filter grows a deeper torso. At `M100`, Q does the opposite violent move: Stage 2 leaps from 980 Hz to 4.75 kHz while Stage 6 falls from 933 Hz to 412 Hz, stretching the anatomy top and bottom at once. The harshest transition is `M0_Q0 -> M100_Q0`: peak level jumps about 26 dB, dynamic range expands about 25 dB, and the centroid drops about 1.27 kHz as the bright shell suddenly acquires a resonant mouth and chest. It feels physical because one fixed 3 kHz face stays put while the jaw, chest, and roof exchange pressure around it.

### What Would Break It
Stage 6 is the load-bearing piece. Without that 225-1103 Hz floor, the filter loses its chest and turns into a bright comb mask with no body under the 3 kHz faceplate.

### Character Tags
[metallic] [vocal] [pressure] [comb] [morphing]

### Relatives
Closest to P2k_013 because both are mislabeled phasers that actually speak through chest, mouth, and bite zones. P2k_012 is less vowel-specific and more comb-scarred: it keeps a fixed 3 kHz mask while other stages cross-swap around it. It also shares some upper-shell behavior with P2k_000, but P2k_000 stays more fractured; P2k_012 holds together as a body with a metal face.
