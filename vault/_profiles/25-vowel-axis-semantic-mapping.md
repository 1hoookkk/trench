# 25-Vowel Axis: Semantic Mapping (The Q-Axis Rule)

## 1. Axiom: Q is Not Just Resonance
In conventional DSP, Q (Quality factor) controls the bandwidth of a specific filter stage. In TRENCH doctrine—specifically when reverse-engineering the physical modeling behavior of E-mu's complex Z-plane structures—**Q is a structural allocator.** The Q-axis controls the literal custody and spatial arrangement of physiological resonance bands (Sub, Chest, Throat, Mouth, Presence, Bite, Air).

It does not simply sharpen a peak; it completely shifts the semantic meaning of the filter structure across its quadrant.

## 2. The Custody Trade (From Autopsies like P2k_016)
When navigating the `M0_Q0` -> `M100_Q100` 4-corner bilinear interpolation space, the Morph axis typically drives the primary vowel transition (e.g., "Ah" -> "Ee"). 

However, the Q axis dictates **how the stages divide the labor to achieve that vowel**. As Q transitions from 0 to 100%, stages systematically trade structural custody.

**Example Matrix Behavior (Low Morph `M0` at `Q=0` vs `Q=100`):**
*   `M0_Q0`: Stages 2 and 3 sit bundled together near 500 Hz (Throat), while Stage 6 pulls up to 1.5 kHz. The "mouth" is physically narrow, creating a specific vowel color.
*   `M0_Q100`: The structure violently realigns. Stage 6 plummets to 357 Hz (Anchoring at Sub/Chest). Stage 2 jumps to 1.3 kHz (Mouth). Stage 3 leaps to 2.3 kHz (Presence).

**The Rule:** Raising Q spreads the formants anatomically, fracturing single dense bands into widely spaced anchor+shingle arrays.

## 3. The Broken Quadrant Phenomenon
Because Q redistributes structural roles, the behavior at High Morph + High Q (`M100_Q100`) often departs entirely from vocal physical modeling and crosses into "fracture" territory.

*   At `M100_Q0`, the stages still adhere to a physical vocal tract sequence (Chest -> Mouth ladders).
*   At `M100_Q100`, the constraints are intentionally overloaded. Floor anchors (Chest/Throat) are preserved, but upper stages (Mouth/Presence) are violently sheared upward into 10–14 kHz "roof" bands. 

## 4. Implementation Constraints for Authoring
When authoring new bodies for the 25-Vowel array or any semantic mapping set, you must obey the following:
1.  **Anchor Preservation:** Do not move the Stage 6 anchor out of the Sub/Chest/Throat band when modulating Q. If the floor is lost, the psychoacoustic illusion of a "vocal cavity" dissolves into a disjointed comb filter.
2.  **Corner Purity:** Because TRENCH strictly enforces 4-corner bilinear interpolation without mid-points, the "vowel sweep" (Morph axis) and the "structural shatter" (Q axis) rely purely on the absolute endpoints. Therefore, you must author `Q0` for maximum phonetic clarity and `Q100` for maximal physical spread/stress.
3.  **No Teleportation:** Formants should smoothly traverse or split from neighboring bands (e.g., Throat -> Air, Mouth -> Presence). Do not map a sub-band directly into a presence band without bridging, or the interpolation will sweep through a dead zone.

## 5. Artifact Manifestation
This behavior maps directly into the runtime `compiled-v1` cartridge:
*   `M0_Q0` arrays contain closely grouped biquad `[c0..c4]` coefficients targeting 400-2000 Hz.
*   `M100_Q100` arrays contain coefficients representing steep notches and starkly separated peaks, leveraging the full 20 Hz – 18 kHz spectrum to simulate material/tract failure.

**Conclusion:** The semantic mapping of a TRENCH filter body requires treating the X-axis (Morph) as the desired target state and the Y-axis (Q) as the mechanical stress/band distribution method.
