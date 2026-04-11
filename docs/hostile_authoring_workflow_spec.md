# TRENCH Hostile Authoring Workflow Spec

Status: proposed main workflow
Owner: product / authoring
Scope: TRENCHWORK body authoring UX and solve loop

## 1. Purpose

Replace direct filter editing with a dead-simple, hostile, absurdist ritual workflow that is still serious enough to author lush, expressive bodies for TRENCH.

The user should not think in:
- corners
- stages
- zeros
- radius
- frequency bytes

The user should think in:
- what must survive
- what must happen across the sweep
- what kind of physical failure the body should enact

This workflow becomes the default authoring surface.

## 2. Core Product Position

TRENCHWORK is not a coefficient editor for normal people.
It is a ritual machine that turns silly acts into serious body trajectories.

The absurdity is front-end language.
The seriousness is:
- bounded latent intent
- deterministic target generation
- hard runtime gates
- solver-backed export

This workflow has a dual mandate:
- be the best authoring method for iconic TRENCH bodies
- be naturally legible as short-form content and product marketing

If it succeeds only as a joke, it fails.
If it succeeds only as an internal tool, it fails.

## 3. Dominant UX Read

The user is completing stupid little hostile tasks to coerce a sound into existence.

The experience should read as:
- punitive
- funny
- tactile
- fast
- musically consequential

It must not read as:
- a synth editor
- a node graph
- a scientific plotting tool
- a game UI
- a generic creativity app

## 4. Product Bar

This workflow is only valid if it beats the current direct-edit experience on both axes below.

### 4.1 Authoring Bar

The workflow must make it easier to author bodies with:
- stronger invariants
- more legible landmark motion
- faster keep/reject iteration
- less solver misuse
- less accidental HF takeover
- less false confidence from pretty plots

Success is not:
- users having fun clicking things
- a cool demo video
- a funny prompt list

Success is:
- better bodies in less time
- fewer invalid candidates
- faster convergence toward iconic presets

### 4.2 Marketing Bar

The workflow must also produce:
- immediately readable on-camera interaction
- memorable prompts that are easy to quote
- audible before/after moments within seconds
- a distinct TRENCH identity that competitors would not invent

The content angle is not decoration.
It is a distribution advantage built into the authoring surface.

## 5. Strategic Rule

Every ritual must satisfy both of these:
- it maps to a real sonic verb
- it is funny or striking enough to be watchable in a short clip

If a ritual is only funny, delete it.
If a ritual is only technically meaningful, hide it in debug.

## 6. Main Workflow

The main workflow is a short ritual sequence.

1. Choose a body brief
2. Complete 5-7 ritual prompts
3. Hear the solved candidate update after each step
4. Press `KEEP`, `REJECT`, or `AGAIN`
5. Save/export only after runtime gates pass

Target completion time:
- one pass in under 20 seconds
- one meaningful keep/reject decision in under 30 seconds

## 7. Ritual Surface

Each screen shows:
- one giant imperative sentence
- one crude visual token or icon
- one primary interaction area
- live audio preview
- one small state caption

Each screen must allow only one of these interaction types:
- repeated tapping
- press-and-hold
- drag through resistance
- choose one of three cursed options

No default sliders.
No raw numbers.
No visible DSP terms.

## 8. Example Ritual Vocabulary

These prompts are not jokes pasted on top of knobs. Each one maps to a stable perceptual verb.

Examples:
- `Pick 10 apples`
- `Paint the door`
- `Solve the captcha`
- `Kick the vending machine`
- `Tear the cardboard slowly`
- `Kiss the dog`
- `Drop the floor out`

Verb families:
- organize / solve / count = coherence, articulation, readability
- paint / smear / coat = breadth, coverage, body, darkness
- kick / tear / break = instability, choke, fracture, collapse
- hold / carry / protect = invariant preservation

## 9. Compute Architecture

The workflow is:

`ritual input -> gesture features -> latent intents -> sweep checkpoints -> global solve -> runtime gates -> audition`

The user never edits solver parameters directly.

### 9.1 Gesture Features

Rituals produce measurable features such as:
- tap count
- tap spacing
- hold duration
- drag distance
- drag speed
- drag wobble
- reversals
- overshoot
- release snap
- hesitation time

These are normalized to `0..1`.

### 9.2 Latent Intent Vector

Gesture features map into a bounded intent vector.

Minimum intent schema:
- `mass`
- `choke`
- `rip`
- `rattle`
- `cry`
- `fracture`
- `coherence`
- `stability_budget`
- `aggression`

Optional later additions:
- `darkness`
- `spread`
- `wetness`
- `betrayal`

Rules:
- intents are bounded
- intents are perceptual, not mathematical
- rituals may affect multiple intents
- each ritual should have one dominant intent and one secondary intent

### 9.3 Sweep Checkpoint Targets

The intent vector populates explicit target checkpoints across the morph sweep.

Default checkpoint grid:
- `m=0.0`
- `m=0.2`
- `m=0.4`
- `m=0.6`
- `m=0.8`
- `m=1.0`

Each checkpoint stores target descriptors such as:
- dominant peak bands
- desired notch bands
- low-band anchor requirement
- top-end clamp or release
- comb density
- allowed instability
- local coherence

Q behavior is not an implied resonance sweep.
Q behavior is owned by the world-level gate contract and only biased by the body route.

Reference:
- `docs/world_q_behavior_gate_spec.md`

### 9.4 Solver

The solver owns the contradiction.

Inputs:
- a chosen pole skeleton or seed body
- checkpoint targets
- invariant constraints
- runtime safety limits

Outputs:
- four runtime corners
- keyframe JSON
- compiled JSON
- diagnostic fit report

The solver may operate in the TRENCH-native four-corner domain.
The user should not see corners during authoring.

## 10. Speaker Knockerz Mapping

The first body this workflow should support is `Speaker Knockerz`.

Invariant:
- sub below ~60 Hz remains anchored across the sweep
- body remains bass-centric
- no HF takeover
- no midpoint or open-state blowups

Default target landmarks:
- `m=0.0`: Vault
- `m=0.2`: Chest Resonance
- `m=0.4`: Choke
- `m=0.6`: Cardboard Rip
- `m=0.75`: Rattle
- `m=0.9`: Cone Cry
- `m=1.0`: Total Fracture

Suggested ritual mapping:
- `Pick 10 apples` -> mass, darkness, floor anchor
- `Paint the door` -> chest coverage, spectral coating, low-mid spread
- `Kick the vending machine` -> choke pressure, instability onset
- `Tear the cardboard slowly` -> rip placement, rip sharpness, raggedness
- `Solve the captcha` -> rattle articulation, coherence, local peak legibility
- `Kiss the dog` -> cone cry focus and late-sweep vocal precision
- `Drop the floor out` -> fracture severity within safety limits

## 11. Example Ritual-to-Intent Mappings

### 11.1 `Paint the door`

Meaning:
- coat a surface
- increase body coverage
- smear edges

Primary gesture features:
- stroke count
- total coverage
- speed
- missed spots

Primary intents:
- `mass`
- `darkness`
- `spread`

Checkpoint effect:
- heavier early sweep
- wider low-mid body
- less surgical notching

### 11.2 `Solve the captcha`

Meaning:
- force order out of disorder
- make an event readable

Primary gesture features:
- completion accuracy
- retries
- completion time

Primary intents:
- `coherence`
- `cry`
- `rattle`

Checkpoint effect:
- narrower, more articulate late landmarks
- less diffuse breakup

### 11.3 `Kiss the dog`

Meaning:
- localized tenderness under absurd pressure

Primary gesture features:
- hold duration
- release snap
- wobble

Primary intents:
- `cry`
- `aggression`
- `stability_budget`

Checkpoint effect:
- late narrow scream
- bounded local violence

## 12. UI Layers

### 12.1 Primary Layer

Visible by default:
- ritual instruction
- interaction zone
- current audio preview
- keep/reject control
- small state phrase

### 12.2 Secondary Layer

Visible on demand:
- checkpoint names
- intent bars
- gate state

State labels should be verbal, not numeric:
- `HEAVIER`
- `MORE RAGGED`
- `TOO CLEAN`
- `THROAT HOLDING`
- `RATTLE FORMING`

### 12.3 Debug Layer

Hidden behind explicit action:
- `SHOW THE MACHINE`

Debug reveals:
- checkpoint targets
- current latent intents
- runtime corner plots
- solver residuals
- gate failures
- exact export payload

The debug layer is for power use only.
It must not pollute the main surface.

## 13. Render Doctrine

The visual reference is Source Engine / Garry's Mod prop-jank.

The target is:
- flat lighting
- hard shadows
- low-poly props
- cheap concrete rooms
- checkerboard or bargain-basement industrial materials
- ugly red utility HUD text
- obvious asset reuse
- slight scale wrongness

The rooms should look like:
- someone dragged 2006-era props into a concrete test box
- nothing expensive
- nothing premium
- nothing sleek
- nothing intentionally “beautiful”

The style is:
- uncanny cheapness
- accidental-looking badness
- deliberate readability

### 13.1 Material Rule

Use:
- grey concrete
- stained metal
- flat painted surfaces
- cheap plywood
- crude warning colors

Do not use:
- premium metals
- soft luxury materials
- glossy sci-fi panels
- ornamental texture detail

### 13.2 Lighting Rule

Use:
- flat overhead lighting
- hard-edged shadows
- minimal bounce
- no cinematic grading

Do not use:
- bloom
- volumetrics
- dramatic fog
- painterly color grading
- moody hero lighting

### 13.3 Asset Rule

Use:
- chunky low-poly objects
- simple collision silhouettes
- visible repetition
- slightly wrong proportions
- obvious room kits

Do not use:
- bespoke hero props in every room
- highly detailed meshes
- realistic clutter simulations
- decorative set dressing with no interaction meaning

### 13.4 Composition Rule

Each room must have:
- one dominant prop
- one obvious interaction zone
- dead space around the prop
- one readable failure state

The prop must dominate the read immediately.

### 13.5 Anti-Rule

Do not build:
- a nostalgia recreation
- a branded Source clone
- a meme map
- a polished “retro-inspired” look

The point is not homage.
The point is cheapness as pressure.

## 14. Safety and Gate Model

The workflow must be hostile in language but safe in compute.

Hard gates before export:
- midpoint audit
- open-state peak not HF-dominated
- sub anchor preservation
- no absurd gain peaks
- world-level Q behavior contract

If a ritual drives the hidden target outside a safe region:
- clamp the target
- preserve the joke
- do not expose the clamp
- annotate internally that the solver hit a safety boundary

The user should feel resistance, not see a validation form.

Q failures should be surfaced internally and in debug as body-language failures, not math jargon.

Examples:
- `Q MAX COLLAPSED TO ONE PEAK`
- `Q LOST THE BODY`
- `Q PUSHED THE SOUND INTO HF`

## 15. Keep / Reject Loop

The sift loop is mandatory.

Controls:
- `KEEP`
- `REJECT`
- `AGAIN`

Behavior:
- `KEEP` saves candidate and exact target state
- `REJECT` discards candidate and resamples around the current latent intent
- `AGAIN` keeps the ritual state but requests another solve nearby

This is the main taste loop.
The solver proposes.
The user judges.

## 16. Short-Form Content Constraint

This workflow must work on camera.

Requirements:
- every ritual is readable in less than 1 second
- every ritual causes an audible change
- every screen has one obvious action
- the absurd prompt should be funny even without DSP context
- the result should still feel like serious sound design

The product pitch should be naturally demonstrable as:
- “I made this body by painting a door and solving a captcha”

## 17. Kill Criteria

Kill the workflow, or redesign it, if any of these happen:
- the average kept body is not better than the current direct workflow
- users start opening debug immediately because the rituals are too vague
- the prompts are memorable but produce muddy sonic results
- the prompts are technically valid but visually dead on camera
- the ritual list becomes random meme sludge with no stable sonic grammar

## 18. Non-Goals

Do not make this:
- a corner editor
- a modular patcher
- a continuous spline authoring tool
- a response-curve drawing app
- a fake game with decorative scoring

Do not expose:
- `a1`
- `r`
- `val1`
- `val2`
- `val3`
- raw corner coordinates
- raw stage selectors

Those belong in debug or engineering tools only.

## 19. Relationship To Heritage

This workflow is not heritage-parity authoring.

Heritage parity remains a separate workflow for:
- RE validation
- compiler parity
- corpus regression

This workflow is TRENCH-native.
It may use heritage bodies or cubes as seeds, but the main interface is perceptual and ritualized.

## 20. RE Translation Layer

The workflow should translate clean-room RE truth into authoring constraints, seeds, and solver priors.
It must not translate heritage UI or historical product metaphors into the front-end ritual surface.

### 20.1 What To Translate From RE

Use RE-derived data for:
- runtime surface truth
- solver constraints
- seed skeletons
- validation gates
- body priors

Concrete structural truths to translate:
- four-corner runtime basis
- Q-first bilinear interpolation
- stage-count and bypass behavior where historically relevant
- historical sample-rate families where parity work matters
- corpus-learned stage activity patterns
- known pole/zero vocab families from sanitized corpus data
- known failure modes such as midpoint blowup and open-state drift

RE truth belongs in:
- solver constraints
- hidden checkpoint priors
- debug overlays
- parity validation

It does not belong in:
- ritual text
- room themes
- user-facing interaction labels
- body-world art direction

### 20.2 What Must Stay Original

Do not translate from RE:
- historical UI metaphors
- Morph Designer terminology
- byte-level concepts into the main UX
- heritage naming or branding
- any visible “compiler parity” framing

The user should never feel like they are operating an emulator editor.
They should feel like they are coercing a hostile machine.

### 20.3 Translation Pipeline

The correct one-way flow is:

`clean-room RE truth -> normalized body priors -> body-world target grammar -> solver -> runtime gates -> keep/reject`

This implies:
- RE is structural input
- body worlds are original authoring shells
- solver is the bridge
- runtime gates remain final authority

### 20.4 Body-World Inputs

Each body world should consume a `BodyWorldSource` bundle with:
- invariant definitions
- checkpoint landmarks
- allowed instability ranges
- seed skeleton candidates
- stage-role priors
- gate limits
- optional corpus similarity references

The room layout and rituals are then authored on top of that bundle.

### 20.5 Example: Speaker Knockerz

For `Speaker Knockerz`, the RE translation layer should provide:
- allowable seed skeletons derived from validated current-runtime bodies or safe heritage scaffolds
- proven sub-anchor-preserving stage patterns
- known dangerous open-state and midpoint regimes
- landmark-friendly frequency regions for:
  - Vault
  - Chest Resonance
  - Choke
  - Cardboard Rip
  - Rattle
  - Cone Cry
  - Total Fracture

The body world then translates those into rooms and rituals.

Example:
- RE says a certain scaffold holds sub and survives interpolation
- the world turns that into `Hold The Piston`
- the solver uses the scaffold as its starting skeleton

### 20.6 Shared Engine, Per-World Grammar

The system architecture should be:
- one shared corridor engine
- one shared solver bridge
- one shared gate layer
- one shared debug layer
- one world per invariant / failure family
- one body route per authored sentence inside that world

This is the correct split.

Do not build one generic world for all bodies.
Do not build separate engines per body.
Do not build one isolated world per flagship body if the bodies share the same gate grammar.

### 20.7 Provenance Rules

Every world should declare:
- which invariants came from product doctrine
- which constraints came from clean-room RE truth
- which seed families came from sanitized corpus analysis
- which room mappings are original TRENCH design inventions

This keeps the boundary clear:
- math may be inherited
- world meaning must be authored

## 21. World / Body / Candidate Architecture

The main workflow should be expressed in three layers:
- one world per invariant / failure family
- one body per authored route inside that world
- one candidate per solver run

Worlds are the grammar.
Bodies are sentences.
Candidates are attempts.

### 21.1 Ownership Split

Worlds own:
- invariant set
- gate vocabulary
- Q behavior contract
- allowed checkpoint types
- latent intent subspace
- seed skeleton families
- failure grammar
- render kit and world material language

Bodies own:
- checkpoint sequence
- which checkpoint types are used
- room order
- room-specific ritual prompts
- route-specific Q bias
- checkpoint weighting
- default seed prior
- audiovisual phrasing and route flavor

Candidates own:
- one concrete ritual run
- one concrete latent intent state
- one concrete checkpoint target set
- one solver result
- one gate result bundle
- one keep/reject outcome

### 21.2 Why This Split

The structural win is that worlds own the gate vocabulary.

Example:
- a `Bass Pressure` world knows what `sub anchor survived` means
- a `Bass Pressure` world knows what `Q max collapsed to one peak` means
- a `Brittle Air` world knows what `mid scoop held` means
- a `Dual Formant` world knows what `two formants remained dominant` means
- a `Null / Fracture` world knows what `state boundary collapse` means

Bodies inside the same world inherit those gates.
They do not need their own duplicated gate stack.

### 21.3 Checkpoint Grammar

Worlds define the allowed checkpoint types.

Bodies define:
- which checkpoint types they use
- in what order
- with what weights

Example:
- `Bass Pressure` world may allow:
  - `anchor`
  - `bloom`
  - `choke`
  - `rip`
  - `rattle`
  - `cry`
  - `fracture`

- `Speaker Knockerz` body may use:
  - `anchor -> bloom -> choke -> rip -> rattle -> cry -> fracture`

- another bass body may use:
  - `anchor -> bloom -> choke -> cry -> fracture`

That keeps the world as a grammar and the body as a sentence.

### 21.4 Shared Engine

Shared across all worlds:
- fixed camera logic
- object state machine framework
- gesture feature extraction
- latent intent accumulation
- checkpoint target builder
- solver bridge
- keep/reject/again loop
- debug layer
- provenance display

### 21.5 Body Route Rule

One room should correspond to one checkpoint family or invariant operation.

The user should not be:
- roaming freely
- solving a puzzle
- managing inventory
- learning lore

They should be:
- entering a room
- enduring one hostile action
- hearing one part of the body change
- moving forward

### 21.6 Camera Rule

Default camera should be:
- fixed
- top-down, shallow isometric, or locked cinematic angle

Do not use free first-person as the main workflow.
First-person adds navigation burden and hurts readability.

### 21.7 Physics Rule

The world may look physically broken.
The authoring mapping must not be physically random.

Use:
- fake strain
- deterministic state machines
- bounded interaction windows
- measured input features

Do not use:
- free sandbox chaos as the source of truth
- uncontrolled rigid-body randomness
- passive failure loops unrelated to sonic intent

### 21.8 Starter World Set

The initial useful world set is:
- `Bass Pressure`
- `Brittle Air`
- `Dual Formant`
- `Null / Fracture`

Each world should launch with one flagship body route first.

### 21.9 Speaker Knockerz Route In Bass Pressure World

`Speaker Knockerz` should be the first full body route.

Recommended room sequence:
- Hold The Piston
- Drag The Crate
- Jam The Gate
- Rip The Sheet
- Sort The Bolts
- Hold The Whistle
- Break The Floor

This route becomes the validation slice for the whole workflow.

## 22. Implementation Phases

### Phase 1: Prototype

Build one-body prototype for `Speaker Knockerz` with:
- one hostile corridor
- 6-7 rooms
- one RE-informed source bundle
- one seed skeleton family
- one solver backend
- keep/reject/again
- hidden debug drawer
- provenance view in debug only

Phase 1 should be reduced further for validation:
- one world: `bass_pressure`
- one body route: `speaker_knockerz`
- first 3 rooms only:
  - `Hold The Piston`
  - `Drag The Crate`
  - `Jam The Gate`

The purpose of this slice is to prove:
- hostile fixed-camera room workflow
- deterministic strain-to-intent mapping
- solver bridge viability
- keep/reject loop speed
- better authoring than direct edit for at least one body family

Implementation brief:
- `C:\Users\hooki\Trench\docs\phase1_bass_pressure_vertical_slice.md`

### Phase 2: Generalize

Add body-specific ritual packs:
- `Bismuth Shrapnel`
- `The Glottal Snare`
- `Phase Molt`

### Phase 3: Adaptive Ritual Sets

Swap ritual wording and task style based on body family while preserving the same intent schema.

## 23. Technical Stack Recommendation

The recommended v1 implementation stack is:
- `@react-three/fiber`
- `@react-three/drei`
- `@use-gesture/react`
- `@react-three/rapier`
- `zustand`
- Web Worker for solver and gate evaluation
- Web Audio API for preview and spatial placement
- Vite for app build and PWA packaging

This stack is recommended for:
- speed of prototyping
- fixed-camera room rendering
- deterministic interaction plumbing
- mobile and browser deployment

### 23.1 Renderer Rule

Do not treat `wgpu/WGSL` as a product requirement for v1.

For v1:
- rendering is disposable
- workflow quality is the actual target
- solver / intent mapping correctness matters more than renderer purity

Use a custom WebGPU / WGSL renderer only if the concept survives and later proves that:
- the look cannot be achieved cleanly in the prototype stack
- the performance ceiling matters
- the render pipeline itself becomes a differentiator

### 23.2 Physics Rule

Physics is presentation and interaction feel.
Physics is not the authoring truth.

Use physics to:
- sell strain
- provide object response
- communicate instability

Do not use physics to:
- decide filter targets directly
- introduce uncontrolled randomness into authoring
- replace deterministic feature extraction

### 23.3 Runtime Boundary

The frontend should be separated from the solver path.

Recommended boundary:
- frontend handles rooms, props, gestures, HUD, and room audio
- worker handles feature normalization, intent mapping, checkpoint generation, solver calls, and gate evaluation

The frontend should speak in:
- `WorldSpec`
- `BodySpec`
- `CandidateSpec`

Not ad hoc local room state only.

## 24. Measurement Plan

The workflow should be measured against the current authoring path.

Authoring metrics:
- median time to first keep
- keeps per 20 candidates
- gate-pass rate before audition
- iconic-body completion rate
- number of debug openings per session

Content metrics:
- hook clarity in first 2 seconds
- average prompt readability on screen
- number of ritual phrases worth clipping
- average audible contrast per ritual step

Translation metrics:
- how often RE-informed seeds outperform unguided seeds
- which body worlds rely most heavily on inherited priors
- whether debug users can trace a kept body back to clear structural inputs

Prototype metrics for the 3-room slice:
- can one pass be completed in under 15 seconds
- does each room create an audible consequence within 1 second
- does `KEEP / REJECT / AGAIN` feel faster than direct parameter nudging
- does the slice produce at least one candidate worth keeping

## 25. Open Technical Questions

- exact latent intent ranges
- checkpoint weighting per body family
- how many retries before solver diversity collapses
- whether `AGAIN` should jitter solver seed, checkpoint weights, or both
- whether one ritual can intentionally consume safety budget from a previous ritual
- how much RE prior should be exposed in debug before it becomes a crutch
- which seed families are stable enough for `Speaker Knockerz` world v1
- whether each body needs its own latent intent subspace or only world-specific mappings

## 26. Acceptance Criteria

- A new user can complete a full authoring pass without seeing raw DSP parameters.
- The default workflow never exposes corners or stage-level editing.
- Each ritual has a stable mapping to bounded latent intents.
- Every candidate is solved through the real runtime path and checked against hard gates.
- `KEEP / REJECT / AGAIN` is fast enough to support rapid audition.
- A debug layer exists for engineering inspection, but the main flow remains absurd and simple.
- The workflow improves authoring outcomes relative to the direct-edit baseline.
- The workflow produces clips or demos that are inherently watchable without extra explanation.
- RE-derived structural truth is consumed by the solver and debug layer without leaking heritage UI into the main workflow.
- The first body world proves the architecture before any second world is built.
