# Cartridge Authoring Schema

## Minimum Required Fields

Every shipping cartridge must specify:

```yaml
name: "string"           # Identity name (not a heritage name)
house: "string"          # Behavioral family
profile: Core | Character | Extreme

body_law:                # What the body IS at rest
  description: "string"  # One sentence: static posture
  invariant_anchor: "string"  # What never changes regardless of morph/Q
  key_frequencies_hz: [f32]   # Where the body lives spectrally
  low_q_identity: "string"    # What you hear at Q=0

behavior_law:            # What the body DOES under modulation
  description: "string"  # One sentence: movement character
  morph_verb: "string"   # What morph does (glides, steps, pivots, scatters)
  q_verb: "string"       # What Q does (tightens, screams, fractures)
  best_driver: "string"  # ENV follower, slow LFO, manual, onset
  coherence_range: [f32, f32]  # Morph coherence bounds

breakage_law:            # How the body FAILS
  description: "string"  # One sentence: failure mode
  break_region: "string" # Where on XY terrain breakage occurs
  break_character: "string"  # What breakage sounds like
  recovery: "string"     # How/whether it comes back

xy_terrain:              # Morph surface behavior
  terrain_family: "string"    # SmoothRamp | CornerVolcano | DiagonalCanyon | FracturedField
  basins: "string"       # Where basins reinforce core truth
  ridges: "string"       # Where ridges expose danger
  dead_zones: "string"   # Known flat/uninteresting regions

dsp_traits:              # Technical fingerprint
  topology: "string"     # Formant | Sculpted | SlopeColor | Focused | Frame | Wah
  stage_count: u8        # Active stages (1-6)
  morph_species: "string"    # Breathe | Tilt | Migrate | Swap | Bloom | Scatter
  peak_dynamic_range_db: f32
  centroid_span_hz: f32

success_criteria:        # When to ship this cartridge
  - "string"             # Concrete testable conditions
rejection_criteria:      # When to kill this cartridge
  - "string"             # Concrete fail conditions
```

## Three Laws

The three laws correspond to the 3-Layer Behavioral Stack (see `00_PRODUCT.md`).

### Body Law — "What body has the signal entered?"
What the body IS. The invariant spectral shape. What you hear at rest with moderate Q and centered morph. This is the cartridge's physics — it doesn't change under modulation, it gets revealed or obscured.

### Behavior Law — "How does that body move under normal use?"
What the body DOES. How it responds to morph sweep, Q intensification, and envelope following. The verb of the cartridge. A body with no behavior law is a static filter preset, not a TRENCH cartridge.

### Breakage Law — "What happens when you break it?"
How the body FAILS. What happens at the terrain's edges — high Q, extreme morph positions, fast modulation. Breakage must be authored, not accidental. If a cartridge can't break interestingly, it may still be a valid Core body, but it cannot be Character or Extreme.

### Naming Rule

The name should evoke the body. The three verbs should explain each layer.
Example: *Pressure Vault = Containment / Ring / Buckle.*
If you can't name the three verbs, the cartridge isn't authored yet.

## Validation Checklist

Before shipping:
- [ ] Body law verified: static presentation is useful at low/mid/high Q
- [ ] Behavior law verified: morph sweep produces audible, coherent change
- [ ] Breakage law verified: extreme positions break in the documented way
- [ ] XY terrain matches documented basins/ridges
- [ ] No verbatim P2K/heritage coefficients (clean-room compliant)
- [ ] Ear-validated on at least 3 source types (vocal, pad, bass)
- [ ] Profile assignment confirmed by listening, not just metrics
