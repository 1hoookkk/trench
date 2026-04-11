## Early Rizer (P2k_025, filterType=25)

### What It Is
A violent comb‑mask that behaves like a tube with a hinged skull: it can sit as a dark, hollow “early sweep” full of carved‑out presence voids, or snap into a bright speaking plate with ~20 dB teeth. Morph at low Q drags an entire roof cluster down into the 1.5–5 kHz face while keeping a fixed ~553 Hz throat bone and a fixed ~17 kHz cap; Q can then abruptly turn that cap into a 1.5 kHz mouth driver. Use it for risers that *talk* (not just open), basses that chew, and transitions where the midrange fractures then re‑seals.

### Stage Roles
- Stage 1: air lid at 4.92–11.82 kHz (sibilance → air band). Both. Stress.
- Stage 2: anchor rung at 0.55–1.95 kHz (F1 → F2 hinge), plus a presence‑scar zero at 2.18–4.16 kHz. Both. Anchor.
- Stage 3: formant mask at 1.51–3.89 kHz (schwa F2 ↔ presence peak), with a ~6.55–7.20 kHz sibilance cut. Both. Identity.
- Stage 4: air roof hinge at 3.87–15.19 kHz, with a fixed ~14.90 kHz air‑notch. Both. Stress.
- Stage 5: bridge plate at 2.47–6.54 kHz (F3 edge ↔ sibilance); its zero either lives at ~16.85 kHz or collapses to ~4.54 kHz to gouge presence. Both. Stress.
- Stage 6: cavity vent at 1.50–16.99 kHz with a near‑unit‑circle notch that jumps ~5.57 kHz ↔ ~17.96 kHz. Both. Cavity.

### Motion Law
This isn’t a smooth “rise”; it’s two coupled machines. With Q low, morph pulls Stages 1/3/4/5 down together (8.8/3.8/11.4/6.2 kHz → 4.9/1.5/3.9/2.5 kHz), packing the face into the telephone/presence zone while Stage 2 stays pinned at 553 Hz and Stage 6 stays pinned at 17 kHz. The result is a pressurized mouth box with brutal hollowing (notches reach ~-138 dB at `M100_Q0`).

Q is the switchblade. At low morph, `M0_Q0 -> M0_Q100` flips Stage 6 from a 17 kHz cap into a 1.50 kHz mouth driver while Stage 2 leaps 0.55 → 1.95 kHz and peaks jump 2.6 → 21.6 dB (centroid 1.76 → 4.62 kHz). It stops being “dark comb” and becomes a bright speaking grate with a deep ~5.57 kHz vent cut.

At high Q, morph becomes a cross‑swap: Stage 2 drops 1.95 → 0.67 kHz as Stage 6 rockets 1.50 → 11.73 kHz and ejects its notch to ~18 kHz (`M0_Q100 -> M100_Q100`: centroid 4.62 → 1.34 kHz, dynamic range +28 dB). That opposing throat‑down / skull‑up motion is the physical illusion: jaw sinks while the headplate lifts.

### What Would Break It
Stage 6. Remove the near‑unit‑circle vent (and its 17 kHz ↔ 1.5 kHz pole flip) and you lose the two‑mode identity; it becomes just a busy high‑mid EQ ladder.

### Character Tags
[vocal] [comb] [violent] [fracture] [morphing]

### Relatives
Closest to `P2k_024` for comb‑jaw intent and to `P2k_020` for the load‑bearing cavity vent; unlike both, `P2k_025`’s signature is the Stage‑6 pole teleport (air cap → mouth driver). Opposite of the hooded clamp of `P2k_019`.
