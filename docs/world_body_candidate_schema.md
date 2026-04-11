# TRENCH World / Body / Candidate Schema

Status: proposed
Purpose: define the concrete data model for hostile world authoring
Scope: shared engine, per-world grammar, per-body routes, per-candidate solve state

Related:
- `docs/world_q_behavior_gate_spec.md`

## 1. Hierarchy

The authoring system is expressed in three layers:

- `WorldSpec`
- `BodySpec`
- `CandidateSpec`

Ownership model:
- world = grammar
- body = authored sentence
- candidate = one attempted utterance

## 2. Design Rules

### 2.1 WorldSpec owns

- invariant set
- gate vocabulary
- Q behavior contract
- allowed checkpoint types
- latent intent axes
- seed skeleton families
- failure grammar
- render kit

### 2.2 BodySpec owns

- checkpoint order
- checkpoint weights
- room sequence
- ritual prompts
- route-specific Q bias
- default seed priors
- route-specific audiovisual flavor

### 2.3 CandidateSpec owns

- one interaction run
- extracted gesture features
- latent intent values
- concrete checkpoint targets
- solver request and result
- gate result bundle
- keep/reject outcome

## 3. Schema

## 3.1 WorldSpec

```ts
type WorldId =
  | "bass_pressure"
  | "brittle_air"
  | "dual_formant"
  | "null_fracture";

type CheckpointType =
  | "anchor"
  | "bloom"
  | "choke"
  | "rip"
  | "rattle"
  | "cry"
  | "fracture"
  | "sheen"
  | "fold"
  | "void"
  | "formant_a"
  | "formant_b"
  | "boundary";

type IntentAxis =
  | "mass"
  | "darkness"
  | "spread"
  | "coherence"
  | "stability_budget"
  | "aggression"
  | "choke"
  | "rip"
  | "rattle"
  | "cry"
  | "fracture"
  | "wetness"
  | "betrayal";

type GateKey =
  | "midpoint_audit"
  | "open_state_balance"
  | "sub_anchor"
  | "gain_ceiling"
  | "q_profile"
  | "mid_scoop"
  | "dual_formant_count"
  | "boundary_collapse";

interface InvariantRule {
  id: string;
  description: string;
  measurement_key: string;
  hard_limit: number | string;
}

interface GateDefinition {
  key: GateKey;
  description: string;
  fail_condition: string;
  metric_keys: string[];
}

interface QStageExpectation {
  phase: "q_min" | "q_mid" | "q_max";
  description: string;
}

interface QHardGate {
  metric_key: string;
  comparator: "<=" | ">=" | "=" | "!=";
  value: number | string | boolean;
  fail_label: string;
}

interface QBehaviorProfile {
  summary: string;
  intent: string[];
  allowed_behaviors: QStageExpectation[];
  forbidden_failures: string[];
  hard_gates: QHardGate[];
}

interface SeedFamilyRef {
  id: string;
  description: string;
  source_kind: "heritage" | "trench" | "corpus_cluster";
  source_refs: string[];
}

interface RenderKit {
  camera_mode: "fixed_topdown" | "fixed_isometric" | "fixed_cinematic";
  material_theme: string;
  room_palette: string[];
  object_style: string;
  hud_style: string;
  lighting_style: "flat_overhead_hard_shadow";
  cheapness_rules: string[];
  banned_visual_traits: string[];
}

interface WorldSpec {
  id: WorldId;
  name: string;
  description: string;
  invariants: InvariantRule[];
  gates: GateDefinition[];
  q_behavior: QBehaviorProfile;
  allowed_checkpoint_types: CheckpointType[];
  latent_intents: IntentAxis[];
  seed_families: SeedFamilyRef[];
  failure_grammar: string[];
  render_kit: RenderKit;
  provenance: {
    re_inputs: string[];
    doctrine_inputs: string[];
    notes: string[];
  };
}
```

## 3.2 BodySpec

```ts
type BodyId = string;

interface CheckpointUse {
  id: string;
  type: CheckpointType;
  order: number;
  morph_pos: number;
  q_bias: number;
  weight: number;
  description: string;
}

type RoomVerb = "hold" | "drag" | "mash" | "sort";

interface RitualPrompt {
  room_id: string;
  prompt_text: string;
  verb: RoomVerb;
  dominant_checkpoint_id: string;
  dominant_intent: IntentAxis;
  secondary_intents: IntentAxis[];
}

interface RoomSpec {
  id: string;
  label: string;
  checkpoint_id: string;
  object_name: string;
  verb: RoomVerb;
  failure_state_label: string;
  state_phrase_pool: string[];
}

interface BodySpec {
  id: BodyId;
  world_id: WorldId;
  name: string;
  invariant_overrides?: Partial<Record<string, number | string>>;
  q_bias_notes?: string[];
  checkpoint_route: CheckpointUse[];
  rooms: RoomSpec[];
  ritual_prompts: RitualPrompt[];
  default_seed_family_ids: string[];
  route_notes: string[];
  marketing_hook: string;
}
```

## 3.3 CandidateSpec

```ts
type CandidateVerdict = "pending" | "keep" | "reject";

type GestureFeatureKey =
  | "tap_count"
  | "tap_spacing_mean"
  | "tap_spacing_variance"
  | "hold_duration"
  | "drag_distance"
  | "drag_speed_mean"
  | "drag_reversals"
  | "overshoot"
  | "correction_count"
  | "wobble"
  | "release_snap"
  | "hesitation_time";

interface GestureFeatureMap {
  [key: string]: number;
}

interface IntentState {
  [key: string]: number;
}

interface CheckpointTarget {
  checkpoint_id: string;
  descriptors: {
    peak_bands?: Array<{ hz: number; weight: number; width: number }>;
    notch_bands?: Array<{ hz: number; depth: number; width: number }>;
    comb_density?: number;
    instability?: number;
    coherence?: number;
    low_anchor?: number;
    top_clamp?: number;
  };
}

interface SolverRequest {
  world_id: WorldId;
  body_id: BodyId;
  seed_family_id: string;
  intent_state: IntentState;
  checkpoint_targets: CheckpointTarget[];
}

interface SolverResultRef {
  keyframe_json_path?: string;
  compiled_json_path?: string;
  diagnostics_path?: string;
}

interface GateResult {
  passed: boolean;
  failures: string[];
  measurements: Record<string, number | string>;
}

interface CandidateSpec {
  id: string;
  world_id: WorldId;
  body_id: BodyId;
  run_index: number;
  room_feature_log: Record<string, GestureFeatureMap>;
  intent_state: IntentState;
  checkpoint_targets: CheckpointTarget[];
  solver_request: SolverRequest;
  solver_result: SolverResultRef;
  gate_result: GateResult;
  verdict: CandidateVerdict;
  notes?: string;
}
```

## 4. World-Level Grammar

Worlds define the allowed checkpoint vocabulary.
Bodies choose a route through that vocabulary.

Example:

`bass_pressure` allows:
- `anchor`
- `bloom`
- `choke`
- `rip`
- `rattle`
- `cry`
- `fracture`

`Speaker Knockerz` uses:
- `anchor -> bloom -> choke -> rip -> rattle -> cry -> fracture`

A second bass body might use:
- `anchor -> bloom -> choke -> cry -> fracture`

This keeps the world as grammar and the body as sentence.

## 5. Example WorldSpec

```json
{
  "id": "bass_pressure",
  "name": "Bass Pressure",
  "description": "Sub-anchored stress bodies that squeeze low mass into choke, rip, rattle, and cry states.",
  "invariants": [
    {
      "id": "sub_anchor",
      "description": "Sub below ~60 Hz remains anchored across the sweep.",
      "measurement_key": "worst_sub_drop_db",
      "hard_limit": 6.0
    },
    {
      "id": "no_hf_takeover",
      "description": "Open-state peak must remain bass-centered or at least not upper-HF dominated.",
      "measurement_key": "open_state_peak_band",
      "hard_limit": "not_hf_dominant"
    }
  ],
  "gates": [
    {
      "key": "midpoint_audit",
      "description": "Interior morph points must stay below the approved gain ceiling.",
      "fail_condition": "worst_peak_db > 35.0",
      "metric_keys": ["worst_peak_db", "worst_point"]
    },
    {
      "key": "sub_anchor",
      "description": "Sub anchor must survive the sweep.",
      "fail_condition": "worst_sub_drop_db > 6.0",
      "metric_keys": ["worst_sub_drop_db", "worst_sub_point"]
    },
    {
      "key": "gain_ceiling",
      "description": "No absurd absolute gain spikes.",
      "fail_condition": "worst_peak_db > 35.0",
      "metric_keys": ["worst_peak_db"]
    },
    {
      "key": "q_profile",
      "description": "Q max must stay bass-bodied and must not collapse to one isolated peak.",
      "fail_condition": "dominant_peak_count_qmax < 2 or stage_domination_ratio_qmax > 0.70 or hf_takeover_flag_qmax = true",
      "metric_keys": [
        "dominant_peak_count_qmax",
        "occupied_bands_qmax",
        "stage_domination_ratio_qmax",
        "hf_takeover_flag_qmax"
      ]
    }
  ],
  "q_behavior": {
    "summary": "Q tightens pressure and sharpens choke, rip, rattle, and cry without collapsing the body to one whistle.",
    "intent": [
      "preserve sub authority",
      "increase articulation",
      "permit pinch without identity loss"
    ],
    "allowed_behaviors": [
      {
        "phase": "q_min",
        "description": "Broad and heavy, bass-led."
      },
      {
        "phase": "q_mid",
        "description": "Clearer choke and rip articulation while still body-led."
      },
      {
        "phase": "q_max",
        "description": "Aggressive and pinched, but at least two occupied regions remain."
      }
    ],
    "forbidden_failures": [
      "single narrow upper-mid or HF laser peak",
      "DC pile satisfying all low landmarks",
      "one-stage domination of the audible read"
    ],
    "hard_gates": [
      {
        "metric_key": "dominant_peak_count_qmax",
        "comparator": ">=",
        "value": 2,
        "fail_label": "Q MAX COLLAPSED TO ONE PEAK"
      },
      {
        "metric_key": "occupied_bands_qmax",
        "comparator": ">=",
        "value": 2,
        "fail_label": "Q LOST THE BODY"
      },
      {
        "metric_key": "stage_domination_ratio_qmax",
        "comparator": "<=",
        "value": 0.7,
        "fail_label": "ONE STAGE TOOK OVER Q MAX"
      },
      {
        "metric_key": "hf_takeover_flag_qmax",
        "comparator": "=",
        "value": false,
        "fail_label": "Q PUSHED THE BODY INTO HF"
      }
    ]
  },
  "allowed_checkpoint_types": [
    "anchor",
    "bloom",
    "choke",
    "rip",
    "rattle",
    "cry",
    "fracture"
  ],
  "latent_intents": [
    "mass",
    "darkness",
    "spread",
    "coherence",
    "stability_budget",
    "aggression",
    "choke",
    "rip",
    "rattle",
    "cry",
    "fracture"
  ],
  "seed_families": [
    {
      "id": "p2k_006_family",
      "description": "Current-runtime-safe bass pressure scaffold family derived from P2k_006 behavior.",
      "source_kind": "heritage",
      "source_refs": ["cartridges/p2k/P2k_006.json"]
    }
  ],
  "failure_grammar": [
    "pressure build",
    "fundamental separation",
    "paper rip",
    "late scream",
    "controlled collapse"
  ],
  "render_kit": {
    "camera_mode": "fixed_isometric",
    "material_theme": "scuffed_concrete_bass_chamber",
    "room_palette": ["#3A3A3A", "#707070", "#991100", "#111111"],
    "object_style": "chunky_industrial",
    "hud_style": "ugly_red_utility",
    "lighting_style": "flat_overhead_hard_shadow",
    "cheapness_rules": [
      "checkerboard or bargain-basement industrial flooring",
      "nothing expensive or premium",
      "obvious asset reuse",
      "slight scale wrongness",
      "one dominant prop per room"
    ],
    "banned_visual_traits": [
      "bloom",
      "volumetrics",
      "sleek sci-fi panels",
      "premium metals",
      "nostalgia cosplay",
      "branded Source clone look"
    ]
  },
  "provenance": {
    "re_inputs": [
      "4-corner bilinear runtime",
      "current-runtime-safe seed family",
      "shipping invariant gates"
    ],
    "doctrine_inputs": [
      "SHIPPING.md Speaker Knockerz brief",
      "TRENCH identity as stress simulator"
    ],
    "notes": [
      "World grammar is original TRENCH design."
    ]
  }
}
```

## 6. Example BodySpec

```json
{
  "id": "speaker_knockerz",
  "world_id": "bass_pressure",
  "name": "Speaker Knockerz",
  "q_bias_notes": [
    "Q should tighten choke, rip, rattle, and cry.",
    "Q must preserve bass body and at least two occupied regions at q_max."
  ],
  "checkpoint_route": [
    {
      "id": "vault",
      "type": "anchor",
      "order": 0,
      "morph_pos": 0.0,
      "q_bias": 0.0,
      "weight": 1.0,
      "description": "Clean sub, heavy weight, clamped highs."
    },
    {
      "id": "chest",
      "type": "bloom",
      "order": 1,
      "morph_pos": 0.2,
      "q_bias": 0.1,
      "weight": 0.8,
      "description": "Warm 150 Hz chest bloom."
    },
    {
      "id": "choke",
      "type": "choke",
      "order": 2,
      "morph_pos": 0.4,
      "q_bias": 0.3,
      "weight": 1.0,
      "description": "Pressure notch near 200 Hz."
    },
    {
      "id": "rip",
      "type": "rip",
      "order": 3,
      "morph_pos": 0.6,
      "q_bias": 0.2,
      "weight": 1.0,
      "description": "Peaky unstable resonance near 400 Hz."
    },
    {
      "id": "rattle",
      "type": "rattle",
      "order": 4,
      "morph_pos": 0.75,
      "q_bias": 0.4,
      "weight": 0.9,
      "description": "Comb-like structural flutter."
    },
    {
      "id": "cry",
      "type": "cry",
      "order": 5,
      "morph_pos": 0.9,
      "q_bias": 0.6,
      "weight": 1.0,
      "description": "Localized late scream near 1.2 kHz."
    },
    {
      "id": "fracture",
      "type": "fracture",
      "order": 6,
      "morph_pos": 1.0,
      "q_bias": 0.7,
      "weight": 0.8,
      "description": "Mid collapse with sub preserved."
    }
  ],
  "rooms": [
    {
      "id": "hold_piston",
      "label": "Hold The Piston",
      "checkpoint_id": "vault",
      "object_name": "hydraulic_piston",
      "verb": "hold",
      "failure_state_label": "FLOOR SLIPPING",
      "state_phrase_pool": ["HEAVIER", "FLOOR HOLDING", "TOO LOOSE"]
    },
    {
      "id": "drag_crate",
      "label": "Drag The Crate",
      "checkpoint_id": "chest",
      "object_name": "concrete_crate",
      "verb": "drag",
      "failure_state_label": "BODY THINNING",
      "state_phrase_pool": ["CHEST BUILDING", "TOO CLEAN", "MORE WOOD"]
    }
  ],
  "ritual_prompts": [
    {
      "room_id": "hold_piston",
      "prompt_text": "HOLD THE FLOOR DOWN",
      "verb": "hold",
      "dominant_checkpoint_id": "vault",
      "dominant_intent": "stability_budget",
      "secondary_intents": ["mass", "darkness"]
    },
    {
      "room_id": "drag_crate",
      "prompt_text": "DRAG THE WEIGHT HOME",
      "verb": "drag",
      "dominant_checkpoint_id": "chest",
      "dominant_intent": "mass",
      "secondary_intents": ["spread", "coherence"]
    }
  ],
  "default_seed_family_ids": ["p2k_006_family"],
  "route_notes": [
    "Bass-centric across the sweep.",
    "No HF takeover."
  ],
  "marketing_hook": "I built this 808 destroyer by holding down a piston and tearing a wall."
}
```

## 7. Example CandidateSpec

```json
{
  "id": "speaker_knockerz_run_0007",
  "world_id": "bass_pressure",
  "body_id": "speaker_knockerz",
  "run_index": 7,
  "room_feature_log": {
    "hold_piston": {
      "hold_duration": 0.84,
      "release_snap": 0.11,
      "wobble": 0.09
    },
    "drag_crate": {
      "drag_distance": 0.73,
      "drag_speed_mean": 0.46,
      "drag_reversals": 0.12,
      "overshoot": 0.18
    }
  },
  "intent_state": {
    "mass": 0.72,
    "darkness": 0.61,
    "spread": 0.39,
    "coherence": 0.48,
    "stability_budget": 0.82,
    "aggression": 0.31,
    "choke": 0.66,
    "rip": 0.58,
    "rattle": 0.52,
    "cry": 0.44,
    "fracture": 0.37
  },
  "checkpoint_targets": [
    {
      "checkpoint_id": "vault",
      "descriptors": {
        "low_anchor": 0.89,
        "top_clamp": 0.72
      }
    },
    {
      "checkpoint_id": "chest",
      "descriptors": {
        "peak_bands": [
          {"hz": 150.0, "weight": 0.68, "width": 0.41}
        ],
        "coherence": 0.51
      }
    }
  ],
  "solver_request": {
    "world_id": "bass_pressure",
    "body_id": "speaker_knockerz",
    "seed_family_id": "p2k_006_family",
    "intent_state": {
      "mass": 0.72,
      "darkness": 0.61
    },
    "checkpoint_targets": []
  },
  "solver_result": {
    "keyframe_json_path": "cartridges/candidates/Speaker_Knockerz_v7.keyframe.json",
    "compiled_json_path": "cartridges/candidates/Speaker_Knockerz_v7.compiled.json",
    "diagnostics_path": "artifacts/speaker_knockerz_run_0007.json"
  },
  "gate_result": {
    "passed": true,
    "failures": [],
    "measurements": {
      "worst_peak_db": 30.2,
      "worst_sub_drop_db": 2.8
    }
  },
  "verdict": "keep",
  "notes": "Strong chest and clean cry without HF takeover."
}
```

## 8. Engine Implications

The shared corridor engine should be driven by these schemas:

- loading a `WorldSpec` defines the available grammar and gate layer
- loading a `BodySpec` defines the route and ritual wording
- each run emits a `CandidateSpec`

This avoids:
- duplicated gate logic
- duplicated invariant logic
- duplicated seed-family logic
- duplicated room engine code

## 9. Prototype Stack Boundary

Recommended v1 frontend stack:
- `@react-three/fiber`
- `@react-three/drei`
- `@use-gesture/react`
- `@react-three/rapier`
- `zustand`
- Web Worker for solver and gate path
- Web Audio API for preview

Rules:
- renderer is disposable in v1
- schema and workflow are not disposable
- physics is presentation and feature extraction support, not authoring truth

Runtime split:
- frontend owns room state, object presentation, gesture capture, HUD, and audio placement
- worker owns intent mapping, checkpoint targets, solver calls, gate evaluation, and candidate persistence

## 10. First Implementation Target

The first implementation should be:
- `WorldSpec`: `bass_pressure`
- `BodySpec`: `speaker_knockerz`
- `CandidateSpec`: persisted for every run
- first 3 rooms only:
  - `Hold The Piston`
  - `Drag The Crate`
  - `Jam The Gate`

Do not build a second world until the first one proves:
- better kept bodies
- better keep/reject speed
- better on-camera readability
