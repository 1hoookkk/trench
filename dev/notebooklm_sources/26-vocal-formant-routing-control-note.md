# Vocal Formant Routing Control Note

Source type: interpretive routing note

Status: usable as a control and labeling framework for vowel-style filters

Constraint: do not treat this as a validated global engine law

## Core Claim

For the dedicated vocal formant family, the control routing can be understood as a two-axis system:

- Axis 1: vowel sweep
- Axis 2: physical-shape metaphor

This is useful as a patching and labeling framework.

It is not yet strong enough to be promoted into a universal doctrine about how the second axis behaves for every body.

## Two-Axis Routing Mechanic

### Axis 1: Vowel Sweep

Routing a modulation source to `FilFreq` moves the filter through vowel states.

Example description:

- `Ah` toward `Ee`
- `Oo` toward `Ah`

In this reading, the horizontal motion is not best understood as a generic cutoff sweep. It is a vowel-state sweep.

### Axis 2: Physical Metaphor

Routing a second modulation source to `FilRes` or `Q` acts like a physical-shape control metaphor for the vocal family.

Useful semantic interpretations:

- body size
- cavity size
- mouth
- throat opening
- vocal tract scale

The point is perceptual control, not generic synth-style resonance language.

## Product Rule

Use this mapping for vocal-style filters.

Do not overlearn it into:

- `Q always means mouth cavity size`

That leap is unsupported.

The correct split is:

- vocal bodies: second axis can be described as body size or cavity scale
- other bodies: second axis should be described by whatever perceptual stress law actually fits that body

## Body-Specific Second-Axis Meanings

Possible mappings by body family:

- vocal bodies: body size, cavity, mouth, size
- acid bodies: pressure, bite
- bass bodies: weight, collapse
- carve or null structures: yield, dominance transfer
- hostile shelves or scrambled machines: whatever perceptual stress law matches the behavior

This keeps the control language aligned with the locked doctrine that the surface should speak in perceptual roles, not raw Q terminology.

## Authoring Guidance

When documenting or authoring a vocal body, describe it as a sentence.

Example:

- Morph moves from Ah to Ee while body size tightens.

This is the correct use of the mapping:

- strong enough for authoring language
- not strong enough for a universal DSP law

## Plugin Labeling Guidance

For vocal presets, do not keep the vertical axis labeled as `Q` if the body behavior is clearly perceptual.

Preferred labels:

- Body
- Cavity
- Mouth
- Size

For non-vocal bodies, relabel the axis to match the actual perceptual job of that body.

## Hard Boundary

Supported as a practical design rule:

- in vowel filters, the second axis behaves more like a physical-shape metaphor than a generic resonance knob

Not supported yet as universal doctrine:

- that the second axis should always be treated as body size across TRENCH or across the original E-MU family

## Safe NotebookLM Interpretation

When using this note in synthesis or product reasoning:

- treat it as a routing and semantic-labeling guide for the vocal family
- keep the provenance warning attached
- avoid restating it as a confirmed engine-wide fact
