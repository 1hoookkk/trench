---
name: feedback-control-semantics
description: User wants body-specific control labels derived from what morph/Q actually do to each body (like Bat-Phaser using "Freq"/"Hz" instead of Morph/Q)
type: feedback
---

Display control names that match what the axis actually DOES to each specific body, not generic "Morph"/"Q".

**Why:** Bat-Phaser (CPhantomBatman) uses "Freq" and "Hz" instead of Morph/Q because the axes map to notch position and notch depth. The control name IS the body's identity.

**How to apply:** Derive labels from scorer metrics — if morph mostly moves centroid, call it "Sweep". If Q mostly controls air shimmer, call it "Shimmer". Store in RuntimeProgram alongside the body. Display in UI.
