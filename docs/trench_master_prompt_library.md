# TRENCH Master Prompt Library: Core Operational Framework

Status: active governing reference
Owner: implementation doctrine
Scope: all agentic work in this repository

## 1. Purpose

This document is the operational authority for agent workflows in TRENCH.

Its job is to:
- keep Domain 1 implementation clean
- stop RE math from leaking into shipping code
- force repo-native verification instead of cargo-cult process text
- preserve context across long authoring or implementation sessions

If another prompt, skill, or handoff conflicts with this file, this file wins unless the root `AGENTS.md` or an explicit user instruction says otherwise.

## 2. Foundational Doctrine: Contextual Sovereignty

Every non-trivial task must declare the active domain before proposing edits.

Conflating domains is a doctrine violation.

### 2.1 Domain 1: TRENCH Implementation

This is the shipping Rust/Python implementation.

Current ground truth:
- topology: 12-stage serial DF2T cascade
- active stages: 6 active + 6 passthrough
- interpolation basis: 4-corner bilinear, Q first then morph
- runtime coefficient payload: `[c0, c1, c2, c3, c4]`
- DF2T mapping: `[b0, b1, b2, a1, a2]`
- runtime form:
  - `y = c0*x + w1`
  - `w1 = c1*x - c3*y + w2`
  - `w2 = c2*x - c4*y`

Do not use Rossum kernel math in Domain 1.

### 2.2 Domain 2: Firmware RE

This is clean-room reverse-engineering truth from hardware, corpus study, and Ghidra work.

RE kernel form may be expressed as:
- `[2 + b1 / b0, 1 - b2 / b0, 2 + a1, 1 - r^2, b0]`

This is not a drop-in representation for Domain 1.

If you apply Domain 2 kernel assumptions directly to Domain 1 coefficient writes, you risk the exact failure class already seen in this repo:
- `c0 = 2.0` blowouts
- false peak alignment
- invalid null behavior

Use Domain 2 only for:
- constraints
- priors
- provenance
- parity validation
- debug overlays

### 2.3 Domain 3: Future Design Choices

This is speculative architecture.

Examples:
- modulation routing not yet shipped
- post-cascade boost or extra dynamics layers
- alternative world authoring shells
- new solver interfaces not promoted into the current stack

Domain 3 must not mutate Domain 1 behavior until explicitly promoted.

## 3. Operational Command Table

The original `npm` command set is wrong for this repo.
TRENCH is Cargo-first at the root.

### 3.1 Mandatory Verification Commands

Use the smallest truthful verification surface first.

| Command | Purpose | Mandatory Frequency |
|---|---|---|
| `cargo test -p <crate>` | crate-local behavior verification | after every implementation or edit in one crate |
| `cargo test --workspace` | cross-crate verification | before promoting changes that span crates or shared contracts |
| `cargo clippy --workspace --all-targets -- -D warnings` | lint and regression gate | after tests pass for code changes |
| `cargo build --workspace` | compile check when shipping/build integrity matters | before release-facing promotion when compilation surface changed |

Null targets remain:
- static `< -120 dB`
- dynamic `< -90 dB`

Do not claim an audible-path or spectral-math change works unless you ran the relevant verification and measured it.

### 3.2 Context Management Commands

Client slash commands are not repo law.
What is repo law is the handoff artifact.

Mandatory policy:
- when session state becomes fragile, write `dev/SESSION_STATE.md`
- after the handoff is written, the operator may run the client's reset command
- do not rely on lossy compaction rituals

Suggested trigger:
- write the handoff once the session becomes long enough that facts are likely to drift
- if the client exposes `/context`, use it as a diagnostic only
- if the client exposes `/clear`, use it only after the handoff is complete

## 4. Vertical Slice Implementation Prompt

Use this template for feature development when the task requires real code.

```xml
<task_instruction>
  <objective>Implement [Feature Name] via Vertical Slice.</objective>

  <domain>1</domain>

  <verification_steps>
    1. Write a failing test in `tests/` or the owning crate's test module.
    2. Run `cargo test -p <crate>` to confirm the failure.
    3. Modify ONLY the necessary files in the owning crate or the directly coupled caller.
    4. Run `cargo test -p <crate>` to confirm the pass.
    5. If the change crosses crate boundaries, run `cargo test --workspace`.
    6. Run `cargo clippy --workspace --all-targets -- -D warnings`.
  </verification_steps>

  <dsp_constraints>
    - Topology: 12-stage serial DF2T cascade.
    - Active set: strictly 6 active + 6 passthrough stages.
    - Math: strictly follow DF2T runtime form.
    - Invariant: coefficients map to [b0, b1, b2, a1, a2].
    - Interpolation: 4-corner bilinear, Q first then morph.
    - Minimalism: no parallel paths, no gain baked into c4, no fudge factors.
  </dsp_constraints>

  <refusal_rails>
    - Do not import Rossum kernel math into Domain 1.
    - Do not invent new abstractions unless the existing code cannot express the change.
    - Do not hide failed fitness or null behavior behind spec rewrites.
  </refusal_rails>

  <next_action>Propose the smallest failing test that proves the vertical slice.</next_action>
</task_instruction>
```

## 5. "Pick Apples" Authoring And Forge Prompt

Use this prompt family when working in authoring, forge generation, or body search.

### 5.1 Strategy Selection

Pick one target body and one strategy:

- `Speaker Knockerz` -> `Cluster Sweep`
- `Ear Bender` -> `HF Wall Fold`
- `Asymmetric Belch` -> `Shelf Tear`
- `12-Pole Laser` -> `Brickwall Collapse`
- `Basin Sculptor` -> `Octave Notch Basin`

### 5.2 Forge Implementation Law

Order of operations:
1. shelf
2. freq
3. peak
4. morph

If you tune `freq` before `shelf`, you are solving the wrong problem.

Zero Coordination Rule:
- identify the two highest-frequency pole anchor stages
- those two stages must cross-target each other
- all remaining stages aim their zeros at the highest pole frequency
- report the alignment frequencies
- reject candidates that fake success through DC pile or non-local maxima

### 5.3 Candidate Generation Prompt

```text
TASK: Generate 10 coefficient candidates.

Domain: 1 for output coefficients, 2 only for priors if needed.
Target Body: [Selected Body]
Applied Strategy: [Selected Strategy]

Vocabulary:
- Use Shelf / Freq / Peak / Morph
- Forbid internal abstractions such as SemanticBody or SlotSpec unless they already exist in code

Output:
- Generate 10 coefficient candidates
- Frame A and Frame B must be intentional gestures
- Report zero-coordination alignment frequencies
- Flag any candidate whose strongest maximum collapses toward DC instead of the intended band
```

## 6. JSON-Structured Design And CMF Prompt

Use this for faceplate or recovered-lineage design generation.

```json
{
  "role": "Industrial Designer, Recovered Heritage Hardware",
  "task": "Generate TRENCH faceplate visual specification",
  "context": "Cartridge-based morph filter. Aesthetic: industrial recovered lineage, not retro-clone.",
  "format": "JSON Block",
  "cmf_specification": {
    "60_30_10_color_rule": {
      "primary_60": "Verdigris copper patina with erosion pitting",
      "secondary_30": "Cyan thumbwheels (#39E0E0) with internal glow",
      "accent_10": "Pink display window"
    },
    "surface_finish": "Dry media blast industrial matte on copper chassis",
    "edge_wear": "Bronze wear at slot bevels",
    "labels": "Bake only FILTER and TYPE text. Do not hallucinate other labels."
  },
  "layout": "5 physical controls: MORPH slider, Q slider, TYPE selector, ENV toggle, REV toggle"
}
```

## 7. Adversarial Audit Configuration

The required skill file lives at:
- `.claude/skills/audit/SKILL.md`

Its job is to attack:
- architecture bloat
- DC pile fraud
- heritage contamination
- optimizer spec laundering
- fake file assumptions

The audit is mandatory before finalizing non-trivial patches, solver changes, or prompt rewrites that could drift doctrine.

## 8. Negative Constraint Refusal Rails

Hard bans:

| Do Not | Use Instead |
|---|---|
| RBJ cookbook character filters | authored gestures inside the 12-stage serial DF2T cascade |
| parallel signal paths | strictly serial cascade stages |
| gain baked into `c4` | correct DF2T feedback math only |
| new dependencies without approval | existing core logic or explicit permission |
| `freq` before `shelf` | define `shelf` first |
| Rossum kernel writes in Domain 1 | maintain Three Domain Separation |
| spec rewrites to hide optimizer failure | report raw failure metrics |

## 9. Document And Clear Handoff Protocol

Use `dev/SESSION_STATE.md` as the authoritative handoff artifact.

Template:

```md
# TRENCH Session Handoff: [Date]

## Current Truths
- [specific implementation decisions made this session]
- Current Domain Focus: [1, 2, or 3]

## Identified Blockers
- [walls hit, failed measurements, or mathematical inconsistencies]

## Next Steps
- [immediate technical requirements for the next session]

## Invariants Status
- [verification actually run and current result]
```

Activation phrase:

```text
Initiate Context Reset. Summarize current Truths, Domain Focus, Blockers, and Next Steps into dev/SESSION_STATE.md. Ensure all current Domain 1 invariants are documented. Once written, the operator may clear the session. Do not use compact-style summarization as a substitute.
```

## 10. Enforcement Summary

If an agent:
- mixes Domain 1 and Domain 2 math
- skips verification
- adds architecture without present need
- rewrites the spec to excuse failed measurements
- or uses fake repo commands

the work is invalid until corrected.
