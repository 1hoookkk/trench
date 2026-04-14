# Houses and Profiles

## House

A House is a behavioral family — a group of cartridges that share spectral ancestry, movement character, or morphological kinship.

Houses are not genre tags. They describe how a family of bodies behaves under modulation, not what music they're "for."

Examples of what might define a House:
- Shared formant topology (vocal/speaking bodies)
- Shared slope character (acid/ramp bodies)
- Shared fracture behavior (contraband/unstable bodies)
- Shared terrain shape (smooth ramp, corner volcano, diagonal canyon)

A House should contain 2–4 cartridges spanning the Core → Character → Extreme range.

## Profile Definitions

Each profile answers one question from the 3-Layer Behavioral Stack.

### Core — "What body has the signal entered?"

The invariant identity that survives any input. What you hear at rest, at moderate Q, at default morph.

- Usable on first touch
- Musical without modulation
- The safest entry point into the House
- Defined by what it reinforces, not what it hides

Core is the basin — the XY terrain's lowest-energy attractor.

### Character — "How does that body move under normal use?"

The everyday verb, the dialect. What emerges when the envelope follower or slow sweep activates the body.

- Identity-leading — it has a name, a verb, a personality
- Best experienced in motion, not at rest
- The reason someone picks this House over another
- Defined by what it reveals under movement

Character is the mid-zone — the terrain where the body speaks.

### Extreme — "What happens when you break it?"

The cult behavior, the failure mode that triggers addiction. What happens when you push past where it was designed to work.

- Dangerous but authored — the breakage is deliberate
- Not randomly unstable — the failure has a shape
- Useful for sound design, resampling, transitional FX
- Defined by how it breaks, not by how it sounds when safe

Extreme is the ridge — the terrain's edge crossings and cliff faces.

## Naming Emerges from the Stack

The three laws generate the cartridge's name. Example:

> **Pressure Vault** = Containment (Core) / Ring (Character) / Buckle (Extreme)

The name should evoke the body. The three verbs should explain what happens at each layer. If you can't name the three verbs, the cartridge isn't authored yet.

## Positioning Rules

1. Every House must have at least one Core cartridge.
2. Character cartridges must sound meaningfully different from Core under modulation.
3. Extreme cartridges must have a recoverable path back to usable territory.
4. A cartridge that is boring at rest AND boring under modulation belongs in no profile.
5. Profile assignment is based on behavior, not subjective intensity.

## Relationship to Body Families (Code)

In `trench-forge/src/scorer/coherence.rs`, `body_family()` classifies generated bodies:

| Code Family | Profile Mapping | Metric Criteria |
|------------|----------------|-----------------|
| Core | Core | `morph_coherence > 0.65 AND q_persistence > 0.55` |
| Character | Character | Neither Core nor Extreme |
| Extreme | Extreme | `morph_coherence < 0.45 OR (terrain_entropy > 0.70 AND morph_coherence < 0.55)` |

These are generation-time heuristics. Final profile assignment is a curatorial decision.

## Terrain Families (Morph Surface Shape)

Terrain family predicts modulation character:

| Terrain | Behavior | Best Profile |
|---------|----------|-------------|
| SmoothRamp | Breathes with dynamics | Core |
| CornerVolcano | Static Q intensification | Core / Character |
| DiagonalCanyon | Chunky/vocal stepping | Character |
| FracturedField | Chaotic under modulation | Extreme |
