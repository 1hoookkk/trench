---
name: Zero trajectory continuity law — proven on 33 original E-mu P2K skins
description: c-domain LERP guarantees analytic continuity of zero positions through fold bifurcations. 594 trajectories tested, zero discontinuities. Interpolation law resolved.
type: project
---

## Law

c-domain LERP guarantees analytic continuity of zero positions through fold bifurcations. This is an empirically validated property of the E-mu byte encoding, not a hypothesis.

## Evidence

**33 P2K skins × 6 stages × 3 morph paths = 594 trajectories.**

- 35 fold bifurcations found (16/33 skins have at least one)
- Zero discontinuities across all 594 trajectories
- Every apparent "step" was a root-labeling artifact — when tracked as a zero pair, positions are continuous through every bifurcation

### Bifurcation anatomy (P2k_016 stage 4, t=0.970→0.971)

```
t=0.970: conjugate pair at 0.658 ± 0.037i  (inter-zero gap = 0.074)
t=0.971: two reals at 0.693 and 0.626      (inter-zero gap = 0.067)
```

Gap shrinks continuously to zero at bifurcation, then grows on other side. Complex plane positions continuous.

### Talking Mouth stage 5 (anti-formant sweep)

```
t=0.00: real zeros -0.73, -0.95 (high-freq cancellation)
t≈0.05: fold — two reals merge into conjugate pair
t=0.50: zeros at +0.08±0.92i (sweeping CCW)
t≈0.95: fold — pair splits back to two reals
t=1.00: real zeros +0.999, +0.994 (near-DC cancellation)
```

## Mechanism

Under c-domain LERP, numerator coefficients are:
- b0(t) = c4(t) — linear in t
- b1(t) = (c0(t)-2) · c4(t) — quadratic in t
- b2(t) = (1-c1(t)) · c4(t) — quadratic in t

The discriminant b1²-4·b0·b2 is a polynomial in t. It crosses zero smoothly (no discontinuity possible for a polynomial). Zero positions, being roots of a quadratic with polynomial coefficients, are continuous functions of t.

## What this closes

1. **Interpolation law: resolved.** c-domain LERP. No Bark warping, no polar decomposition, no parameter-space recalculation. The 10-byte encoding IS the interpolation space.
2. **Catastrophe theory concern: empirically dismissed.** Fold bifurcations exist but the encoding makes them continuous crossings, not cliffs.
3. **Grisso vs parameter-space approach: settled.** E-mu's encoded coefficient approach handles topological bifurcations that parameter interpolation can't represent.

## What this opens

Forge outputs must satisfy the same continuity property. Any body where c-domain LERP produces a discontinuous zero trajectory is invalid by architectural law, not aesthetic judgment. The 33 P2K skins are the validation oracle.

**Why:** Resolves the fundamental question of whether zero placement during morphing is safe. It is — by construction of the encoding.

**How to apply:** Trust c-domain LERP for all interpolation. No special handling needed. Add trajectory continuity as a hard gate in forge validation — any body failing the 594-trajectory test pattern is rejected.
