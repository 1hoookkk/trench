# TRENCH

Filter plugin by 1hook. MORPH / Q / TYPE. nih-plug VST3/CLAP.

## 1. DIRECTIVES (Read Carefully)

- **Direct Execution:** No preamble. No planning documents. Execute the requested change immediately.
- **No Circling:** Focus on one specific deliverable at a time.
- **No Truncation:** You must output fully qualified, copy-pasteable blocks. If you change a function, provide the ENTIRE function. No `// existing code` or `...rest` placeholders.
- **Division of Labor:** I own sound, taste, UX, and product identity. You own code, math, implementation, and verification.

## 2. REPO STRUCTURE

```text
trench-core/    [CRATE] DSP engine (read trench-core/AGENTS.md before editing)
trench-codec/   [CRATE] Minifloat encode/decode
trench-plugin/  [CRATE] nih-plug wrapper, GUI, runtime params
pyruntime/      [PYTHON] Forge, c4 solver, offline authoring
cartridges/     compiled-v1 JSON payloads
tests/          Integration, null-tests against X3 reference audio
tools/          Scripts, batch processors, validation harnesses

3. VERIFICATION & BUILD
Bash

cargo test -p <crate>
cargo clippy --workspace --all-targets -- -D warnings

    Null Targets: Static < -120 dB, Dynamic < -90 dB.

    Prove It: Do not claim a change works unless you have run the code and verified the output. If the audible path or spectral math changes, state exactly what changed and why.

4. THE INVARIANTS (Frozen, Do Not Deviate)

    Topology: 12-stage serial DF2T cascade. 6 active + 6 passthrough.

    Math: Kernel-form [c0..c4] interpolation only.

    Order: 4-corner bilinear interpolation. Evaluate Q-axis first, then morph-axis.

    Runtime: 32-sample blocks, per-sample coefficient ramping.

    Rates: Authoring = 39062.5 Hz. Plugin = 44100 Hz.

5. STRICTLY BANNED

    No Frameworks: Do not invent abstractions not directly derived from proven hardware data.

    No Band-Aids: Do not add compensation layers, arbitrary saturation, or gain baked into c4 to "fix" plots. Fix the underlying math.

    No Parallel Paths: The engine is a rigid serial cascade.

    No RBJ Cookbook: Clean biquads are spectrally flat. Do not use them for character filters.

    No Unapproved Deps: Do not modify Cargo.toml dependencies without asking.

    No Shipping Heritage: P2K extractions are for dev/calibration only, not shipping presets. No E-mu branding.

    No Pole Sanitization: Never sanitize poles without a written, failing test proving instability.

6. SCOPE

**Nothing is off-limits. Measure it. Build it. Ask me if it sounds right.**

7. STOP AND ASK IF:

    The DSP/audible path changes materially.

    A proposed shortcut creates architectural debt.

    You are unsure if a bug is an implementation error or a deliberate "betrayal" built into the engine's design.
