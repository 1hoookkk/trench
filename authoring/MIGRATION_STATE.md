# Migration State

## Current Phase

Phase 1 complete:

- `authoring/` scaffold created
- canonical entry docs created
- canonical doorway named explicitly

Phase 1.1 complete:

- `FILTER_WORK.md` created as the canonical filter-work repo contract
- `C:\Users\hooki\Trench` named as the canonical repo for filter work
- `trenchwork_clean`, `trenchwork`, `trench_re_vault`, and `do-it` treated as
  reference/history trees, not shipping authorities

Phase 1.2 complete:

- old provenance-taxonomy language removed from the active filter-work route
- the active policy is now: TRENCH does not ship E-mu filters
- reference material can inform, test, or reject candidates, but not define a
  shipped cartridge

Phase 2 complete (2026-05-10):

- `forge/` folded into `authoring/`:
  - `forge/compilers/` → `authoring/compilers/`
  - `forge/audition.py` → `authoring/audition.py`
  - `forge/{authoring_studio,boost_normalizer,extreme_carve,forge_cli,
    hand_authored,iconic_author*,mine_variants,morph_grid_scanner,
    phoneme_author_studio,pill,proof_*,render_reference,sift,stage_rca}.py`
    → `authoring/tools/`
  - `forge/{q_laws,diagnostics,scratch}` → `authoring/{q_laws,diagnostics,scratch}`
  - `forge/{HANDOFF_GPT55.md,SEED.md}` → `authoring/_research/`
- `schemas/` folded into `authoring/schemas/`:
  - `trench.authoring_path.cube.v1.schema.json`
  - `trench.compiled.cube_surface.v1.schema.json`
- E-mu clones in former `forge/compilers/` deleted (Peak/Shelf, EMU Designer,
  Morph2, MorphLP/MorphLPX, Z-plane Qlaw, X3 importers, EMU style banks,
  DillusionMan law verifier, summarizers, EmulatorX_utf8.rc).
- `forge/reproduce_heritage.py` deleted (E-mu by name).
- `runtime/trench-core/tests/peak_shelf_authoring.rs` deleted (E-mu authoring
  test in runtime tree).
- `authoring/peak_shelf/` archived to `authoring/_research/peak_shelf/` with
  a research-only README.
- `./check` updated: `tools/parity_null.py` → `authoring/compilers/parity_null.py`.
- `FILTER_WORK.md` updated: "What Belongs Here" no longer splits authoring/forge;
  "Current Active Filter Path" section no longer points at Peak/Shelf.

## Not Done Yet

- root-level loose Python and scratch files (find_814_*, k.py, optimize_*,
  pro_render, scout_heroes, tmp_parse_filters, _ai_context_*, build_error.txt,
  CANDIDATES.txt, handoff*.json, session-ses_*.md, plus duplicates of
  audition.py / phoneme_author_studio.py / forge_phoneme_author_studio.py /
  audit_pack.py / extract_candidates.py / trench_forge_compiler.py) — pending
  user disposition (delete vs move into `authoring/tools/`).
- sonic tables not yet physically unified
- workbench surfaces not yet indexed from one canonical README
- design references not yet consolidated
- historical `do-it` precedent not yet summarized into `history/`
- external tree inventories are not yet written under `authoring/indexes/`
- no automated no-E-mu-shipping release gate exists yet
- `authoring/{heritage,imported_x3,emu_designer,dirty_airlock}` subdirs
  may contain E-mu material that should be archived under `_research/` —
  pending review.

## Next Recommended Steps

1. close out the root-level cleanup once the operator confirms scope
2. add an automated release gate that rejects E-mu/P2K clone candidates
3. write the first TRENCH-native authoring slice (body-goal-driven, not
   heritage-filter-class-driven)
4. update root docs (`CLAUDE.md`, `dev/SESSION_STATE.md`) to point at the
   consolidated structure
5. inventory and triage `authoring/{heritage,imported_x3,emu_designer}`
   per the FILTER_WORK.md veto rules
