# TRENCH World Q-Behavior Gate Spec

Status: proposed
Purpose: define world-level Q behavior as an explicit gate contract
Scope: solver targets, runtime gates, body-route authoring, candidate rejection

## 1. Purpose

This document closes the gap between:
- "Q is the second cube axis"
- "that sounds wrong at Q max"

Without this contract, the solver and the authoring path are free to confuse:
- useful focus
- with single-stage domination
- or total identity collapse

The core rule:
- Q is not a generic resonance knob
- Q is a world-owned behavior grammar
- bodies inherit that grammar and bias it

## 2. Ownership

Ownership split:
- world owns allowed Q behavior families and hard gates
- body owns route-specific Q bias and emphasis
- candidate owns the measured result

Q behavior is not decided ad hoc per solve.

## 3. Definitions

### 3.1 Terms

- `q_min`: low Q edge behavior, evaluated near `Q=0.0`
- `q_mid`: interior Q behavior, evaluated near `Q=0.5`
- `q_max`: high Q edge behavior, evaluated near `Q=1.0`
- `dominant_peak_count`: count of locally prominent peaks above threshold
- `prominent_peak`: local maximum with required prominence over adjacent bands
- `band_occupancy`: number of distinct occupied landmark bands with persistent energy
- `stage_domination_ratio`: strongest single-stage peak contribution divided by total cascade prominence in the same local region
- `hf_takeover`: dominant peak centroid pushed into forbidden HF region for the world
- `single_peak_collapse`: world identity collapsing to one narrow dominant peak when that behavior is not allowed

### 3.2 Required Measurement Discipline

Peak checks must use:
- local maxima
- prominence relative to neighboring bands
- width
- surrounding suppression

Loose "max in band" checks are invalid. They permit DC piles and monotonic slopes to fake the target.

## 4. Generic Q Gate Primitives

These primitives are shared across worlds.

### 4.1 Peak Persistence

Measures whether multiple meaningful peaks survive as Q increases.

Suggested metrics:
- `dominant_peak_count_qmin`
- `dominant_peak_count_qmid`
- `dominant_peak_count_qmax`
- `peak_prominence_db_q*`
- `peak_width_oct_q*`

### 4.2 Band Occupancy

Measures whether the expected regions stay alive.

Suggested metrics:
- `occupied_bands_qmin`
- `occupied_bands_qmid`
- `occupied_bands_qmax`
- `landmark_band_drop_count`

### 4.3 Single-Stage Domination

Measures whether one stage has taken over the body.

Suggested metrics:
- `stage_domination_ratio_qmax`
- `top_stage_share_qmax`
- `top_stage_share_qcurve`

If one stage owns the entire Q-max read, the body has collapsed unless the world explicitly allows that.

### 4.4 Spectral Center Control

Measures whether the center of attention drifts into a forbidden region.

Suggested metrics:
- `dominant_peak_hz_qmin`
- `dominant_peak_hz_qmid`
- `dominant_peak_hz_qmax`
- `spectral_centroid_hz_qmax`
- `hf_takeover_flag_qmax`

### 4.5 Trajectory Continuity

Measures whether Q sharpens the intended events instead of teleporting to a different body.

Suggested metrics:
- `peak_motion_jump_hz`
- `peak_order_swaps`
- `landmark_persistence_ratio`
- `q_curve_discontinuity_score`

## 5. World Profiles

## 5.1 Bass Pressure

Intent:
- Q should tighten pressure, sharpen choke/rip/cry, and preserve sub authority
- Q may reduce breadth
- Q must not collapse the body into one upper-mid or HF laser peak

Allowed Q behavior:
- `q_min`: broad, heavy, bass-led
- `q_mid`: clearer choke/rip articulation, still bass-centric
- `q_max`: pinched and aggressive, but with more than one occupied identity region

Required at `q_max`:
- at least `2` prominent occupied regions
- no upper-HF dominant takeover
- sub anchor still present
- late narrow feature may appear, but cannot erase the body floor

Forbidden failures:
- single narrow peak above the body with everything else dead
- DC or sub pile satisfying all low landmarks alone
- high-Q endpoint becoming a generic whistle
- one-stage domination of the entire audible read

Suggested hard gates:
- `dominant_peak_count_qmax >= 2`
- `occupied_bands_qmax >= 2`
- `stage_domination_ratio_qmax <= 0.70`
- `hf_takeover_flag_qmax = false`
- `sub_anchor_qmax = pass`

## 5.2 Brittle Air

Intent:
- Q should sharpen brittle shear, glass edges, foil scream, and mid-scoop tension
- collapse toward a narrower bright structure is allowed
- complete single-line collapse is still a failure unless the body brief explicitly authorizes it

Allowed Q behavior:
- `q_min`: spread air, brittle sheen, structural width
- `q_mid`: sharper cuts, clearer shards
- `q_max`: tight glass beam with at least one secondary structure still alive

Required at `q_max`:
- one dominant bright structure is acceptable
- at least one secondary supporting region must survive
- mid scoop must remain legible if the body family requires it

Forbidden failures:
- pure sine-like laser with no structural shoulder
- broadband flattening under high Q
- upper-HF spike with no relation to the route landmarks

Suggested hard gates:
- `dominant_peak_count_qmax >= 1`
- `occupied_bands_qmax >= 2`
- `stage_domination_ratio_qmax <= 0.80`
- `mid_scoop_qcurve = pass`

## 5.3 Dual Formant

Intent:
- Q should clarify the vowel, not delete it
- high Q must preserve at least two meaningful formant regions
- the world is invalid if Q max turns the body into one screaming pole

Allowed Q behavior:
- `q_min`: blurred mouth/throat shape
- `q_mid`: clearer formant split and talk legibility
- `q_max`: narrow but still recognizably dual-formant

Required at `q_max`:
- at least `2` prominent peaks
- preferred `3` if the route asks for it
- persistent `F1/F2` spacing within world-specific range
- no single-peak collapse

Forbidden failures:
- one laser peak replacing the vowel
- F1/F2 merging into one line
- Q max moving the body out of the vocal band family

Suggested hard gates:
- `dominant_peak_count_qmax >= 2`
- `dual_formant_spacing_qmax = pass`
- `stage_domination_ratio_qmax <= 0.65`
- `landmark_persistence_ratio >= 0.60`

## 5.4 Null Fracture

Intent:
- Q may intensify discontinuity, voids, combs, and collapse
- more aggressive identity loss is allowed here than in any other world
- even here, collapse must read as fracture, not accidental single-stage takeover

Allowed Q behavior:
- `q_min`: unstable body with visible internal structure
- `q_mid`: deeper voids, sharper boundaries
- `q_max`: controlled fracture, sparse surviving ridges, partial collapse

Required at `q_max`:
- at least one surviving ridge plus one surviving void/boundary signature
- collapse must still exhibit structured failure, not a monotone line

Forbidden failures:
- total flattening
- pure whistle
- random unrelated dominant band with no void/boundary evidence

Suggested hard gates:
- `boundary_collapse = pass`
- `occupied_bands_qmax >= 1`
- `void_signature_qmax = pass`
- `stage_domination_ratio_qmax <= 0.85`

## 6. Body-Level Q Bias

Bodies do not redefine the world contract.
Bodies only bias how strongly they use it.

Body-level controls may specify:
- which checkpoints respond most strongly to Q
- whether Q should emphasize focus, collapse, articulation, or betrayal
- whether Q max is expected to preserve width or sacrifice it

Examples:
- `Speaker Knockerz`: Q should tighten choke, rip, rattle, and cry, but preserve bass body and at least two occupied regions
- a second Bass Pressure body might skip rattle and let Q concentrate choke plus cry instead

## 7. Candidate Rejection Rules

Reject immediately if:
- world-required `dominant_peak_count_qmax` fails
- `stage_domination_ratio_qmax` exceeds world limit
- `hf_takeover_flag_qmax` fails where forbidden
- a route-required band disappears under Q
- Q max produces a body type outside the world grammar

This is a hard reject, not a taste penalty.

## 8. Implementation Contract

The optimizer and the hand-authoring workflow must both consume the same Q contract.

Required integration points:
- world spec stores Q behavior profile
- body spec stores Q bias and per-checkpoint emphasis
- candidate diagnostics report Q-profile measurements
- keep/reject UI can surface plain-language failures such as:
  - `Q MAX COLLAPSED TO ONE PEAK`
  - `Q LOST THE THROAT`
  - `Q PUSHED THE BODY INTO HF`

## 9. Speaker Knockerz Minimum Q Contract

World: `bass_pressure`

At `q_max`, `Speaker Knockerz` must satisfy all of these:
- sub anchor still audible and measured as present
- at least `2` prominent occupied regions remain
- one region may be late and narrow
- the whole body may pinch, but cannot become one isolated laser peak
- no HF takeover
- no single-stage domination above `0.70`

If `q_max` produces one ugly whistle, the candidate is invalid regardless of any other score.

## 10. Practical Consequence

This document does not make the bodies sound good by itself.

It does make one thing explicit:
- "that sounds wrong at Q max" is now a measurable contract failure
- not a debate
- not a renderer issue
- not a prompt issue
