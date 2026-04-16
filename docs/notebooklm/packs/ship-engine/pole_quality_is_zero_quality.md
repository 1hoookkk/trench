---
name: Pole quality IS zero quality
description: Forge scoring metrics (talkingness, ridge, movement) implicitly score zero-derivation quality. No separate zero metric needed. Mixed radii give the zero law both teeth and subtlety.
type: feedback
---

Pole quality determines zero quality under any runtime derivation law.

- Spread poles → derived zeros at meaningful frequencies → articulated spectral fingers
- Clustered poles → clustered zeros → one big cancellation, flat result regardless of zero law
- Tight poles (r>0.98) → sharp peaks → derived zeros create dramatic notches (the law has "teeth")
- Loose poles (r<0.5) → broad humps → derived zeros just flatten slightly (subtlety)
- Mixed radii in one body → dynamic range across the morph surface

**Why:** The user confirmed that forge scoring metrics (talkingness, ridge_prominence, movement) are implicitly zero-quality metrics. No separate zero scorer is needed. This simplifies the generation-curation pipeline.

**How to apply:** When evaluating pole seeds (from Morpheus cubes or forge generation), the existing scorer already captures zero-derivation potential. Focus generation on pole spread + radius variation, not on explicit zero placement. The zero law (whether baked or runtime) inherits quality from the poles.
