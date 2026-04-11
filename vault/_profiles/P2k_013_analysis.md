## Phaser 2 (P2k_013, filterType=13)

`P2k_013` is labeled **"Phaser 2"**, but the data reads like a **mouth/throat resonator** with surgical cavities — not a whooshy swirl. Misclassification is the point here: it *talks*.

### What Makes It Unique
- **A hard "bone" anchor** down at ~160–200 Hz (plus a low-mid ridge ~900–950 Hz). This is the tether that keeps the rest from becoming random comb soup — it's the jaw hinge.
- **Stacked "mouth formants"** clustered in the speech zone: ~1.5 kHz and ~2.1–2.7 kHz. That's schwa-ish throat (≈1500 Hz) + upper mouth/EE territory (≈2300 Hz), with an ear-canal / ring proximity around ≈2700 Hz. This is why it reads as vocal.
- **A deep cavity cut**: notches can hit ~-100 dB at high Q. That's not "EQ"; that's *air being removed from a tube*. The silence between peaks is the throat shape.
- **A "teeth/metal edge" band** around ~4.3–4.8 kHz, plus a high brace up near ~8–10 kHz. That's the bite and the plastic sheen — the part that makes it feel synthetic and tense instead of warm.

### Stage Roles
- Stage 1: **air/brace** at 8.4–10.2 kHz. Moves with both morph+Q. Stress.
- Stage 2: **formant F1** at 195–953 Hz. Massive Q shift (953→195). Identity.
- Stage 3: **formant F2** at 1.5–2.4 kHz. Rises with Q. Identity.
- Stage 4: **bridge** at 2.2–2.7 kHz. Drifts with Q to fill the gap between F2 and bite. Cavity.
- Stage 5: **bite/teeth** at 4.3–4.8 kHz. Tightens with morph. Stress.
- Stage 6: **bone anchor** at 157–1789 Hz. Fixed at low Q, jumps to mid at high Q. Anchor.

### Motion Law
- **Q is the clamp**: it massively increases carve/depth (dynamic range jumps ~+37 to +40 dB across Q), pushing the "center of attention" upward. This is the throat tightening and the mouth narrowing — vowels become more pronounced and more strained.
- **Morph is the jaw/pressure change**: at low Q it mostly feels like **pressure + posture** (not a huge "shift up"), but at high Q morph also **lifts** the vowel energy upward. That's the illusion: one ridge swells while another tightens and the cavity between them deepens.

### What Would Break It
Stage 6 (the bone anchor). Without it, the formants have no floor — they become a floating comb filter with no sense of body or gravity. The anchor is what makes it feel attached to a chest.

### Character Tags
[vocal] [biological] [cavity] [pressure] [misnamed]

### Relatives
Near-identical to Talking Hedz (filterType=36) at the coefficient level — likely the same filter reindexed for the P2K extended bank. Shares formant territory with Vocal Ah-Ay-Ee (P2k_016) and Vocal Oo-Ah (P2k_017), but those are explicit vowel morphs; this one achieves "talking" through cavity physics, not labeled vowel targets.
