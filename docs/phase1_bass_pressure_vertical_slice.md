# Phase 1 Implementation Brief

Status: approved v1 target
World: `bass_pressure`
Body: `speaker_knockerz`
Scope: 3-room vertical slice

Related:
- `C:\Users\hooki\Trench\docs\hostile_authoring_workflow_spec.md`
- `C:\Users\hooki\Trench\docs\world_body_candidate_schema.md`
- `C:\Users\hooki\Trench\docs\world_q_behavior_gate_spec.md`

## 1. Goal

Build the smallest real authoring slice that can prove all of these:
- hostile room workflow is faster than direct stage poking
- strain-to-intent mapping is deterministic
- `KEEP / REJECT / AGAIN` is viable without a long solver stall
- the body grammar can reject bad Q behavior before export

This build is not a toy demo.
It is a product and authoring feasibility test.

## 2. Non-Goals

Do not build:
- all 7 `Speaker Knockerz` rooms
- a generic world editor
- free roaming
- first-person camera
- real physics-driven authoring
- custom WebGPU renderer
- a heritage-parity frontend

Do not expose:
- stage controls
- corners
- `a1 / r / val1 / val2 / val3`
- raw frequencies in the primary flow

## 3. Stack

Frontend:
- `@react-three/fiber`
- `@react-three/drei`
- `@use-gesture/react`
- `@react-three/rapier`
- `zustand`
- Vite

Backend in browser:
- Web Worker for feature normalization, intent mapping, target building, solver calls, gate evaluation
- Web Audio preview bridge

Rules:
- React owns presentation and local interaction state
- Worker owns all authoring truth
- physics is display feel, not truth

## 3.1 Workbench Integration

This slice does not replace the TRENCHWORK shell.

Use it as:
- one launched authoring mode inside the existing workbench
- one dedicated view or window entered from the normal utility chrome

Do not turn the whole product into a 3D corridor.

Correct split:
- workbench shell remains utility-like
- hostile corridor owns the authoring canvas only

## 4. Slice Content

Use only these 3 rooms:
- `Hold The Piston`
- `Drag The Crate`
- `Jam The Gate`

These 3 rooms must be enough to author:
- anchor / mass
- chest / spread
- choke / early instability

This is intentionally incomplete.
The point is to validate the workflow grammar before adding rip/rattle/cry/fracture.

## 5. Camera And Presentation

Camera:
- fixed isometric or fixed cinematic angle
- no free orbit
- no free walk
- optional micro-parallax only

Room presentation:
- one dominant prop
- one interaction zone
- one failure state label
- dead space around the object

Look:
- cheap concrete chamber
- flat overhead light
- hard shadow
- ugly red utility HUD
- low-poly prop-jank

## 6. App Boundary

## 6.1 Frontend Owns

- room rendering
- prop state machine visuals
- gesture capture
- short-term interaction progress
- audio playback controls
- visible state phrases
- `KEEP / REJECT / AGAIN`
- debug drawer toggle

## 6.2 Worker Owns

- feature normalization
- latent intent accumulation
- checkpoint target build
- seed selection
- local solve request
- runtime gate evaluation
- candidate persistence payload
- plain-language failure labels

## 6.3 Solver Boundary

The worker must not start a full global search on every click.

The worker should perform:
- warm-start seed selection from a safe body bank
- local neighborhood mutation
- bounded rescore
- hard gate rejection

Heavy search stays offline.

## 7. File And Module Layout

Suggested frontend app root:
- `trenchwork-hostile/`

Suggested structure:
- `src/app/App.tsx`
- `src/app/store.ts`
- `src/scene/RoomHost.tsx`
- `src/scene/camera.ts`
- `src/rooms/HoldPistonRoom.tsx`
- `src/rooms/DragCrateRoom.tsx`
- `src/rooms/JamGateRoom.tsx`
- `src/objects/Piston.tsx`
- `src/objects/Crate.tsx`
- `src/objects/Gate.tsx`
- `src/ui/Hud.tsx`
- `src/ui/DebugDrawer.tsx`
- `src/audio/preview.ts`
- `src/worker/authoring.worker.ts`
- `src/worker/contracts.ts`
- `src/worker/intent.ts`
- `src/worker/checkpoints.ts`
- `src/worker/again.ts`
- `src/worker/gates.ts`
- `src/worker/seed_bank.ts`

## 8. State Model

The app should hold:
- current `WorldSpec`
- current `BodySpec`
- active room index
- current `CandidateSpec`
- current interaction progress
- preview audio state
- debug drawer state

Use `zustand`.
Do not spread authoring truth across component-local React state.

## 9. Room State Machines

Each room must be deterministic and short.

## 9.1 Hold The Piston

Primary verb:
- `hold`

Visual:
- heavy piston trying to rise through 2 short stress pulses

Measured features:
- `hold_duration`
- `hesitation_time`
- `release_snap`
- `correction_count`

Maps primarily to:
- `stability_budget`
- `mass`
- `darkness`

State machine:
- `idle`
- `press_started`
- `pulse_1`
- `pulse_2`
- `held_success`
- `released_early`
- `complete`

Duration target:
- `2.0 - 2.5s`

Failure label:
- `PISTON LOST`

## 9.2 Drag The Crate

Primary verb:
- `drag`

Visual:
- concrete crate on a short rail with exaggerated resistance

Measured features:
- `drag_distance`
- `drag_speed_mean`
- `drag_reversals`
- `overshoot`
- `wobble`

Maps primarily to:
- `mass`
- `spread`
- `coherence`

State machine:
- `idle`
- `grabbed`
- `under_load`
- `near_target`
- `overshoot`
- `locked`
- `complete`

Duration target:
- `2.0 - 3.0s`

Failure label:
- `BODY THINNING`

## 9.3 Jam The Gate

Primary verb:
- `mash`

Visual:
- rusted sliding gate reopening in 3 short bursts

Measured features:
- `tap_count`
- `tap_spacing_mean`
- `tap_spacing_variance`
- `release_snap`

Maps primarily to:
- `choke`
- `aggression`
- `coherence`

State machine:
- `idle`
- `burst_1`
- `burst_2`
- `burst_3`
- `sealed`
- `missed_window`
- `complete`

Duration target:
- `1.5 - 2.5s`

Failure label:
- `GATE SLIPPING`

## 10. Interaction Rules

Rules for all rooms:
- one dominant interaction only
- one completion path
- one obvious failure label
- no hidden combos
- the room must resolve in under `3s`

Do not add:
- inventory
- scoring
- roaming
- alternate object targets
- permanent fail loops

## 11. Latent Intent Mapping

For the 3-room slice, use this reduced intent set:
- `mass`
- `darkness`
- `spread`
- `coherence`
- `stability_budget`
- `aggression`
- `choke`

All values:
- normalized `0..1`
- clamped
- monotonic under repeatable input

Mapping rule:
- one dominant intent per room
- at most two secondary intents
- no room may move more than 3 intent axes

## 12. Checkpoint Builder

For the vertical slice, build only these checkpoints:
- `vault`
- `chest`
- `choke`

Route:
- `vault` at `m=0.0`
- `chest` at `m=0.2`
- `choke` at `m=0.4`

Checkpoint builder should emit:
- `peak_bands`
- `notch_bands`
- `low_anchor`
- `coherence`
- `instability`
- `top_clamp`

The builder must be deterministic for the same normalized room feature log.

## 13. Seed Strategy

Do not start from random bodies.

Use a seed bank with:
- safe `bass_pressure` scaffolds
- route-specific defaults for `speaker_knockerz`
- precomputed metadata for:
  - sub anchor
  - midpoint behavior
  - open-state balance
  - Q-profile compliance

Seed selection order:
1. body default seeds
2. world safe seeds
3. nearest kept candidate neighborhood

## 14. AGAIN Strategy

`AGAIN` must feel instant enough to preserve flow.

Rules:
- do not rerun unconstrained DE
- do not restart from random
- keep the same room feature log
- keep the same target grammar
- jitter only:
  - seed choice
  - local target weights
  - bounded mutation radius

Suggested behavior:
- `AGAIN #1`: reseed from nearest safe alternative
- `AGAIN #2`: widen local mutation radius slightly
- `AGAIN #3+`: warn in debug that diversity is decaying

Target latency:
- under `1.5s` for visible candidate refresh
- hard ceiling `3.0s`

If the solve exceeds the ceiling:
- reuse previous audio until the new candidate is ready
- never freeze the room UI

## 15. KEEP And REJECT

`KEEP`:
- stores `CandidateSpec`
- stores room feature log
- stores intent state
- stores checkpoint targets
- stores gate result bundle
- stores seed provenance

`REJECT`:
- discards current candidate
- keeps room log in memory only for nearby resample
- does not pollute saved library

## 16. Gate Contract For Phase 1

Minimum hard gates:
- midpoint audit
- sub anchor
- open-state balance
- gain ceiling
- `bass_pressure` Q profile

Q-profile gate must reject:
- `dominant_peak_count_qmax < 2`
- `occupied_bands_qmax < 2`
- `stage_domination_ratio_qmax > 0.70`
- `hf_takeover_flag_qmax = true`

Plain-language failure labels:
- `Q MAX COLLAPSED TO ONE PEAK`
- `Q LOST THE BODY`
- `Q PUSHED THE SOUND INTO HF`
- `SUB FLOOR GAVE OUT`
- `MIDPOINT BLEW UP`

## 17. Debug Drawer

Hidden by default.

Must show:
- world id
- body id
- current room feature log
- intent vector
- checkpoint targets
- chosen seed id
- gate measurements
- gate failures
- candidate verdict

May show later:
- corner plots
- stage plots
- exact payload

Do not expose debug in the main room frame.

## 18. Audio Behavior

The room must always have audible consequence.

Rules:
- each room completion triggers candidate refresh
- preview should update no later than `1s` after room completion
- room object audio should sound colocated with the prop when possible
- no blocking audio graph rebuild on the main thread

## 19. Success Metrics

The slice passes if:
- one full run can be completed in under `15s`
- each room is readable in under `1s`
- candidate refresh after room completion is under `1.5s`
- `AGAIN` feels faster than direct parameter nudging
- at least one kept candidate from the slice survives hard gates and is worth auditioning

The slice fails if:
- solver latency breaks flow
- users open debug before finishing one pass
- rooms feel random
- Q max failure still slips through as a kept candidate

## 20. Immediate Build Order

1. Implement `WorldSpec` / `BodySpec` / `CandidateSpec` loading in the app shell.
2. Implement `Hold The Piston` with feature extraction and worker round-trip.
3. Implement `Drag The Crate`.
4. Implement `Jam The Gate`.
5. Implement reduced checkpoint builder for `vault / chest / choke`.
6. Implement seed bank and local `AGAIN` strategy.
7. Implement phase-1 gates including Q-profile rejection.
8. Add debug drawer.
9. Run timing test.
10. Audition and decide whether the slice is real.

## 21. Cut Rule

If the first 3-room slice does not produce better candidates faster than direct stage poking, stop.

Do not add more rooms, more art direction, or more engine before proving that.
