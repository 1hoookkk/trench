---
name: LLM Compliance Directives
description: Top 10 mandatory rules for all coding agents modifying the TRENCH codebase — serial cascade, zero-forcing, body-skeleton, no Talking Hedz benchmark
type: feedback
---

Mandatory checklist for TRENCH codebase modifications:

1. **NEVER** use "Talking Hedz" as the benchmark for product identity; it is a debug fixture. TRENCH makes original material.
2. **NEVER** implement parallel signal paths; the engine is a rigid **serial cascade** of 6 biquad stages.
3. **NEVER** use "clean" RBJ cookbook biquads for character filters. They are "spectrally flat." Use **Zero-Forcing** for independent zero placement.
4. **NEVER** modify `generate.rs` to generate stages independently. All generation must start from a **Body-Skeleton Plan** with role bounds.
5. **NEVER** add UI overhead, emotional states, or "noise." Instrument-grade restraint — every element earns its place.
6. **NEVER** use the Oscillator SINC kernel table (`0x18065bb70`) for filter damping. The **Filter Q-Law** source is separate and currently unidentified.
7. **NEVER** gate bodies purely on "dramatic movement" or high morph scores. Character can be subtle; movement is a **ranking signal**, not an admission gate.
8. **NEVER** implement 1D (morph-only) interpolation. The runtime must use **4-corner bilinear interpolation** (morph x Q) with the `duplicated_pair` policy.
9. **NEVER** mutate authoring templates in place. Authoring stages compile into **immutable Runtime Programs**.
10. **NEVER** dump raw binary addresses or large code files into context without summarization (Context Efficiency rule).

**Why:** These rules prevent the most common failure modes when AI agents modify DSP codebases. Each rule has a specific incident or architectural reason behind it.

**How to apply:** Check every code change against this list before committing. If a change violates any rule, stop and reconsider.
