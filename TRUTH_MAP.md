# TRENCH — Active Truth Map

One page. No nostalgia. If a path isn't listed under an active section,
it is either archive/non-canonical or queued for deletion.

## Shipping runtime

- **Plugin**: `trench-juce/plugin/` — **sibling git, not in this repo.**
  C++/JUCE 8. Windows standalone + VST3.
- **DSP source of truth**: `trench-juce/plugin/source/TrenchEngine.cpp`
  (referenced by `DOCTRINE.md`, `SPEC.md`). The 12-stage serial DF2T
  cascade and the bilinear cartridge interpolation live there.
- **Rust crates that serve the shipping runtime**: `trench-core` only.
  It exists to codify DSP math (`cascade.rs`), the cartridge schema
  (`cartridge.rs`), and heritage bakes (`hedz_rom.rs`) that the JUCE
  plugin can consume via FFI or a committed C header. If it doesn't
  serve the JUCE plugin, it doesn't belong here.

## Shipping cartridge format

- **Schema**: `cartridge.schema.json` at repo root.
- **Canonical format tag**: `"format": "compiled-v1"`.
- **Keyframe shape**: 4 keyframes (`M0_Q0`, `M100_Q0`, `M0_Q100`,
  `M100_Q100`) × 6 or 12 stages (12 = 6 active + 6 passthrough tail) ×
  `{c0, c1, c2, c3, c4}` direct-form DF2T coefficients + per-keyframe
  `boost`.
- **Rust loader**: `trench-core/src/cartridge.rs::Cartridge::from_json`
  — accepts both the keyframe format above and the array format
  (`corners.M0_Q0 = [[c0..c4], ...]`).
- **Python serializer**: `pyruntime/body.py::Body::to_compiled_json`.

## Verification command

    ./check

Runs, in order: canonical doc set check, cargo type-check the whole
workspace, `pyruntime` import check, null-target test, 41-entry
P2K parity null (hard fail at > -140 dB), Rust cascade parity vs 33
raw skins × 4 corners (hard fail at > -120 dB), and
`cartridge.schema.json` validation across `cartridges/**/*.json`.

## Authoring tools allowed to touch shipping

These write files the Rust crates or the JUCE plugin actually load.
Everything else is a research script.

| Tool | Writes what |
|---|---|
| `tools/extract_emu_filter_params.py` | `vault/_phonemes/heritage_designer_sections.json` — heritage integer grids (read-only data for forge targeting) |
| `tools/bake_hedz_const.py` | `trench-core/src/hedz_rom.rs`, `trench-core/src/hedz_golden.rs` — hardcoded Talking Hedz cartridge + golden impulse vector for the cross-language parity test |
| `pyruntime/forge_shipping.py` + `pyruntime/forge_shipping_finalists.py` | Shipping body candidates under `vault/_shipping_finalists/` |
| `tools/author_sonic_bank.py` | `vault/_shapes/` shape-bank regeneration |
| `tools/parity_null.py` | P2K parity gate consumed by `./check` |

## Archive / non-canonical paths

Kept for history or research; **not in the shipping chain**. None of
these should be imported, linked, or referenced by anything in the
active runtime or the verification command.

- `cartridges/heritage_quarantine/` — compiled biquad floats from
  trenchwork_clean. Known to diverge from the markdown-XML compile
  path (5-stage vs 6-stage runtime semantics per
  `docs/looperator/talking_hedz_x3_surfaces.md`). Do not treat as
  ground truth.
- `dev/notebooklm_sources/` — research narrative md files. Reference
  material for vocal/formant authoring, not runtime data.
- `docs/notebooklm/packs/` — imported NotebookLM bundle packs
  (ground-truth, ship-engine, frontier). Research context; the real
  ground truth is the integer grid in
  `vault/_phonemes/heritage_designer_sections.json`.
- `vault/_atlas/`, `vault/_cracked/`, `vault/_diagnostics/`,
  `vault/_plots/`, `vault/_profiles/`, `vault/_scorecards/`,
  `vault/_triage/`, `vault/_splice_candidates/`,
  `vault/_vocal_pack_2026-04-11/` — authoring scratch and diagnostics.
- `pyruntime/api.py` and everything `forge_*.py` that isn't
  `forge_shipping*` — FastAPI research cockpit. Serves the authoring
  UI, not the shipping runtime.
- `docs/archive/`, `docs/archive_symbol_inventory_2026-04-09.*` —
  snapshot-at-date inventories; historical only.

## Delete next

Claims with evidence. The purge order is safest-first.

| Path | Why it goes | Evidence |
|---|---|---|
| `trench-live/` (whole crate) | Orphan Rust plugin harness. Zero doc mentions. The shipping plugin is `trench-juce/plugin/` (JUCE, C++, sibling repo). Only references in this repo are its own files + workspace `Cargo.toml:2`. | `grep -rln "trench-live" /home/user/trench/` returns three self-references and nothing else. |
| `trench-codec/` | 119-line crate, referenced only by its own `Cargo.toml` and archive inventories. No doc owns it, no other crate depends on it. Verify before deletion with a workspace-wide import check. | Same grep pattern — only self-references plus archive dumps. |
| `xtask/` | 21 lines in `xtask/src/main.rs`. Cargo convention for workspace tooling; check whether anything in `./check` or CI invokes it before deletion. | `grep -rln xtask` shows `BODIES.md` and an archive handoff — probably historical. |
| `trench-core/src/cluster.rs`, `trench-core/src/motor.rs`, `trench-core/src/role.rs`, `trench-core/src/engine.rs` | 800+ lines of Rust that may or may not be consumed by the JUCE plugin. `engine.rs` is exported via `lib.rs:12` but only its own tests use it inside the workspace. `Cluster` / `Motor` / `Role` similarly — audit what the JUCE plugin actually FFI's into before keeping. | Workspace grep + the sibling-repo constraint. |
| `trench-live/src/editor.rs` | **Already deleted** (Phase 3 of the current surgical slice, pending commit). | `ls trench-live/src/` shows only `lib.rs`. |

## Architect Prompts

Questions you can paste to me during TRENCH work. Each one forces me
to do the digging and report back in plain terms — no jargon, no
hedging, no "it depends", no hiding behind complexity. I own the
implementation. You own the decision. These prompts put me in the
position where I have to give you something you can actually act on.

**What actually runs.** Tell me every file the shipping plugin
actually uses when it makes sound. For each one, show me the
evidence — what loads it, what calls into it. Anything you can't
prove is loaded is either dead or lying. Give me a short list, not
a tour.

**Which path ships.** For `<thing>`, walk me through exactly what
happens from the moment I touch it to the moment sound comes out.
If there's more than one way the codebase could do it, name every
way and tell me which one is actually the one that ships. You are
not allowed to answer "keep both".

**What's actually enforced.** Pick a rule or invariant the code
claims in a comment, a docstring, or a name. Tell me whether
there's a real check — a test or a runtime guard — that enforces
it. If it only exists as words, mark it as a wish, not a rule.
Wishes are liabilities.

**What just broke.** Something failed. Before you write a single
sentence of explanation, do these five things in order:
1. What changed — every commit since the last time the repo was
   green, and which files each one touched.
2. Which code the broken thing actually uses — name the path from
   the plugin entry point to where the failure is.
3. What the failure was — the exact name of the test or check, what
   it expected, what it got, and how far off.
4. Which files are now suspect — the overlap between "files that
   changed" and "files the broken thing uses".
5. What to do — for each suspect, give me one word: REVERT,
   ISOLATE, DELETE, or KEEP, plus one sentence of why.
Only after those five are answered are you allowed to tell me the
root cause story.

**What's safe to delete.** Find every file, folder, crate, or big
function in this repo that nothing else actually uses. For each
one, mark it DEAD (delete it), QUARANTINE (move to an archive
folder so history survives), or LIVE (prove it with one caller).
Sort by size, biggest first. No "maybe useful later" — that's how
the repo got bloated in the first place.

**Where the rules get broken.** The repo has a rule: `<quote it>`.
Check every piece of code that could cross that rule. For each
violation, show me the file and a short description, and tell me
whether it's a bug to fix or whether the rule is fiction. Pick per
violation. No softening, no averaging.

**Where the repo is bluffing.** Name every place this codebase is
claiming more solidity than it has actually earned. I want examples
like: tests that don't really test what their name says, invariants
that live in comments but never get checked at runtime, docs that
describe what we wish the system did as if it already works that
way, claims of correctness that lean on "trust me". One line per
bluff. Rank by how badly it would hurt if I believed it.

**Where the rules are fiction.** Read every rule in `CLAUDE.md`,
`DOCTRINE.md`, `pyruntime/CLAUDE.md`, and `trench-core/CLAUDE.md`.
For each rule, count how many times the actual code already breaks
it. If a rule is broken five or more times, it's not a rule — it's
decoration. Tell me per rule: rewrite the rule to match reality, or
fix the code to match the rule.

**What `./check` actually proves.** Go through every step in
`./check` and tell me, in one sentence each, what the step
genuinely proves — compared to what its name makes it sound like it
proves. Flag every step whose name oversells what the test is
actually locking down.

**What could glitch the audio.** Look at the audio callback and
everything it calls. Find every place in there that allocates
memory, takes a lock, touches a file, parses text, or does anything
that's not pure math. Any of those can make the audio stutter.
Target is zero. For each one, show me where it is and tell me
whether it has to stay or can be killed.

**What could crash.** Find every place in the shipping code that
can blow up when something unexpected happens. For each one, tell
me: is this on purpose (a real invariant, documented, keep it), a
missing sanity check (we should handle this properly), or lazy
(someone skipped error handling and left a bomb)? Show me where.

**Where the same thing lives twice.** Find every concept in the
repo that has more than one implementation. For each: name the
concept, name every place it lives, tell me which one should be
the canonical version, and show me how the copies have drifted
apart. If none of them is obviously canonical, flag it as a
conflict and recommend which one I should make canonical.

**Can this ship.** Given the current state of the repo, list
everything standing between here and a release I can be confident
in. For each one, tell me whether it's a real bug, a test we're
missing, a claim we haven't proven, or architectural debt. Rank by
how bad it would be if we shipped anyway. Nothing "nice to have" —
if it's not blocking the release, it doesn't belong on this list.

**Force a choice.** `<some feature>` has more than one
implementation. List them. Show me the evidence of which one is
live, which is dead, and which is contested. Then tell me: keep
exactly one, delete the rest, one sentence of reasoning per
verdict. You are not allowed to say "keep both", "it depends", or
"it's complicated".

**What this branch shouldn't contain.** This branch was supposed to
do `<X>`. List every change on the branch that isn't `<X>`. For
each stray change tell me: does it belong (fold in), should it move
to a separate branch (split out), or is it a mistake (revert)?
Group them into those three buckets so I can decide the whole pile
at once.

## Notes

- The **sibling-repo invariant** makes everything in this repo's
  `trench-live/`, `trench-codec/`, `xtask/`, and most of `trench-core/`
  suspect until someone audits what `trench-juce/plugin/` actually
  FFI's into. That audit has to happen before the deletion list above
  can be executed safely.
- `trench-core/src/cartridge.rs` + `cascade.rs` + `agc.rs` +
  `hedz_rom.rs` + `hedz_golden.rs` are the only Rust files I can
  defend as shipping-relevant without seeing the sibling repo.
