## Morph Designer (P2k_022, filterType=22)

### What It Is
This is a “design” filter in the physical sense: a bolted-on faceplate with a hinged roof and a vented throat, capable of turning a sound into a serrated mask full of near-silence gaps. At `M100_Q0` it’s brutally carved (peak ~22.3 dB, notches ~-123.1 dB, ~75 peaks/75 notches): a comb-grill you can hear air pressure fighting through. Raise Q and it doesn’t just “resonate more” — it *re-wires* where the roof lives and how sealed the cavity is, swinging from sibilant teeth to a tighter speaking band.

### Stage Roles
- Stage 1: bite/air lid at 2.83–7.75 kHz (presence ↔ sibilance). Both. Stress.
- Stage 2: formant/mask post at 2.96–3.24 kHz (F3/telephone-high zone). Stays put. Identity.
- Stage 3: mask rail at 2.21–3.09 kHz (ear-canal ↔ F3 edge). Both. Identity.
- Stage 4: formant carrier at 1.35–2.54 kHz (F2 mouth). Both. Identity.
- Stage 5: cavity vent at 0.99–1.91 kHz, with a zero that can pin DC (`0 Hz`) or jump up to ~4.89 kHz. Both. Cavity.
- Stage 6: anchor at 227–435 Hz (chest ↔ F1 floor), backed by an air-band zero (~7.58–17.96 kHz). Both. Anchor.

### Motion Law
The anatomy is a stack with one nailed “face” (Stage 2 ~3 kHz) and a roof (Stage 1) that decides whether you’re hearing sibilant metal or a presence mask. At Q0, morph is a jaw-and-crown counter-rotation: Stage 1 climbs 3.48 → 5.11 kHz while Stages 3–5 sink (2.54/1.92/1.56 kHz → 2.21/1.35/0.99 kHz) and the low anchor rises (227 → 435 Hz). That drop in centroid is huge (~-1233 Hz) and the carve gets more violent (range +26 dB).

Q is the hinge release. At `M0`, raising Q lifts *almost everything* upward (Stage 1 launches 3.48 → 7.75 kHz) and partially fills the voids (range -27.5 dB as notches rise from ~-112 dB to ~-75.9 dB). At `M100`, raising Q flips the direction of the roof: Stage 1 dives 5.11 → 2.83 kHz while Stages 3–5 climb together into a tighter 1.9–3.1 kHz mask and the anchor drops (435 → 243 Hz). The most violent transition is `M100_Q0 -> M100_Q100`: peak drops ~13.4 dB and dynamic range collapses ~40 dB as the grill vents and the whole thing snaps from “teeth over a hollow box” into a more coherent talking plate.

### What Would Break It
Stage 5. If you weaken that vent/cavity stage (the one that can literally put a zero at DC), the filter loses the near-silence gaps and reads like a busy midrange EQ instead of a pressurized comb-mask.

### Character Tags
[vocal] [comb] [surgical] [pressure] [morphing]

### Relatives
Closest in *motion logic* to `P2k_012` (fixed ~3 kHz faceplate with cross-swaps around it) and to `P2k_020`–`P2k_021` for “comb voids that assemble into speech.” Opposite of the lid-style lowpass bodies (`P2k_001`–`P2k_003`), which read as a roof closing rather than a vented mask reconfiguring.

