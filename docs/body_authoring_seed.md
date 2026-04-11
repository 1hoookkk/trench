# TRENCH Body Authoring Seed Prompt

Use this prompt in a fresh session when generating or editing the four shipping bodies.

---

## Context

You are authoring 4 shipping filter bodies for TRENCH, a cartridge-based morph filter instrument. Each body is a 12-stage serial DF2T cascade at runtime with:
- `6` active stages
- `6` passthrough stages
- `4` corners
- kernel-form coefficient interpolation only

The runtime is frozen. The authoring layer may use friendlier abstractions, but it must compile back into the existing 4-corner, 6-stage `compiled-v1` format.

Read these first:
- [SPEC.md](C:\Users\hooki\trench-juce\SPEC.md)
- [SHIPPING.md](C:\Users\hooki\trench-juce\SHIPPING.md)

## The 4 Bodies

Only author these names:

| # | Name | Role |
|---|------|------|
| 1 | Speaker Knockerz | Sub weapon. Pressure, choke, cone cry. |
| 2 | Aluminum Siding | Brittle top-end damage. Sibilant fold, foil tear. |
| 3 | Small Talk | Vocal cavity morph. Open throat to strangled bite. |
| 4 | Cul-De-Sac | Structural fracture body. Tube to comb collapse. |

Do not invent alternate public names. Match the filter to the name.

## Method: Intent First, Spectrum Proves It

Do not start by asking how many stages to use or how dense the filter should be. Complex filters are easy. The task is authored intent.

For each body, define:
1. `Identity`
2. `Invariant`
3. `Danger Zone`
4. `Motion Law`

Then prove those claims on the spectrum and only then compile into stage structure.

### 1. Identity

What is the body trying to do on real material?

Examples:
- `Speaker Knockerz`: squeeze bass air mass without killing the floor
- `Aluminum Siding`: weaponize brittle treble without filling the mids
- `Small Talk`: force an anatomical vocal read out of non-vocal material
- `Cul-De-Sac`: cross a structural boundary and return as a different filter

### 2. Invariant

What must survive the full morph?

Examples:
- pinned sub below `60 Hz`
- permanent `~1 kHz` acoustic void
- exactly two dominant formants
- root-note hum tether

If the invariant disappears, reject the candidate.

### 3. Danger Zone

Where does the body become special?

Examples:
- `200 Hz` choke
- `400 Hz` cone rip
- `5 kHz / 8 kHz` sibilant fold
- total null boundary

These are not accidents. They are the authored reason the body is worth using.

### 4. Motion Law

What does movement mean?

Define:
- what `MORPH` exposes
- what `Q` intensifies or betrays
- what, if any, envelope-follow behavior should eventually matter

If motion does not reveal new behavior, the body is under-authored.

## Sonic Tables

Use the sonic tables as guide rails, not decoration.

They should drive:
- formant targets
- anti-resonance targets
- instrument/body landmarks
- notch territories
- stress regions

Use the tables to support intent:
- vowels for cavity identity
- nasals for zero placement
- physical landmarks for material read
- HF landmarks for brittle damage and boundary behavior

Do not use the tables as a preset list. Use them as acoustic vocabulary.

## Structural Method

1. Build or choose an initial spectral scaffold.
2. Place poles and zeros to support the written identity and invariant.
3. Verify the path across all four corners.
4. Only then reduce the result into the shipping 6-stage active runtime.

Hidden authoring complexity is allowed if it improves the final result. Shipping runtime complexity is not.

That means you may:
- explore more internal resonant events during generation
- use helper structures while searching
- generate multiple candidates per name

But the final shipping body must still compile to:
- `compiled-v1`
- `4` corners
- `6` active stages per corner

## Candidate Standard

Generate four candidates per shipping name. Keep only candidates that respect the written brief in [SHIPPING.md](C:\Users\hooki\trench-juce\SHIPPING.md).

Reject a candidate if:
- it loses the invariant
- the danger zone is weak or missing
- the morph path is static or ornamental
- it only works as a frozen sweet spot
- it sounds complex without sounding intentional

## Verification Loop

For each candidate:
1. Inspect the spectrum at all four corners.
2. Check whether the invariant survives.
3. Check whether the danger zone is audible and placed correctly.
4. Check whether morph and Q produce meaningful movement.
5. Compare against the shipping name brief, not just the pretty plot.

If trench MCP is available, use it for:
- stability scans
- state comparison
- corridor / identity drift checks
- ranking the four candidates per family

## Critical Rule

The user is not shipping "a complicated filter."
The user is shipping:
- a bass choke weapon
- a brittle air destroyer
- a talking throat
- a structural fracture body

Author toward those ends directly.
