# Script Inventory

_Generated 2026-04-17_


## C:/Users/hooki/Trench

### Python (22 scripts)

| Path | Purpose | Last Modified | Referenced By | Writes Files |
|------|---------|---------------|---------------|--------------|
| `C:/Users/hooki/Trench/dev/sift/gen_inventory.py` | Build script_inventory.md for Trench, trench_re_vault, trenchwork_clean. | 2026-04-17 | no | yes — C:/Users/hooki/Trench/dev/sift/script_inventory.md |
| `C:/Users/hooki/Trench/tools/bake_hedz_const.py` | Bake Talking Hedz heritage template → Rust const arrays + golden impulse response. | 2026-04-14 | no | yes — heritage_designer_sections.json |
| `C:/Users/hooki/Trench/tools/compile_grid.py` | Compile an authored 2-token grid into a compiled-v1 cartridge. | 2026-04-16 | no | yes — token_inventory_unified_v2.json, {key}.json |
| `C:/Users/hooki/Trench/tools/compile_raw.py` | Compile an internal raw stage surface into a compiled-v1 cartridge. | 2026-04-17 | yes — Trench\tools\forge_write.py | yes |
| `C:/Users/hooki/Trench/tools/extract_emu_filter_params.py` | Extract E-mu MorphDesigner filter parameters from NotebookLM markdown exports. | 2026-04-14 | no | yes |
| `C:/Users/hooki/Trench/tools/forge_write.py` | Forge write wrapper — compile raw-stage-v1 and atomically land it in the | 2026-04-17 | no | yes — forge-watch/active.json |
| `C:/Users/hooki/Trench/tools/run_shipping_release_gate.py` | Fail-closed shipping release gate for emotional-vocal confidence. | 2026-04-13 | yes — Trench\xtask\src\main.rs | yes — better_than_emu_holdout_scorecard.json, shipping_release_gate.json |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/_graveyard/forge_trajectory.py` | Trajectory-first forge — author stress narratives, solve for corners. | 2026-04-07 | no | yes — {key}_{i:03d}.json, trajectory_sweep.png |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/calibration_to_body.py` | calibration_to_body.py — Convert calibration JSON directly to compiled-v1 body. | 2026-04-07 | no | yes |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_generator.py` | forge_generator.py — Anchor + Vector body generator for TRENCH Sift workflow. | 2026-04-07 | yes — Trench\trench-juce\forge\pyruntime\forge_shipping.py | yes |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_heritage_optimize.py` | Heritage-space pymoo optimizer v2 — morph-only search, pressure-derived Q. | 2026-04-07 | no | yes — taste_model.json, {name}.compiled.json |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_optimize.py` | Forge optimizer — pymoo multi-objective body search. | 2026-04-07 | yes — Trench\trench-juce\forge\pyruntime\forge_heritage_optimize.py | yes — taste_model.json, {name}.json |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_plot.py` | Plot frequency responses for TRENCH bodies. | 2026-04-07 | no | no |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_shipping.py` | Forge shipping — P2K direct-synthesis optimizer for the 4 flagship bodies. | 2026-04-07 | no | yes — {name}.json, best_morph_sweep.png |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/recipes/author_weekend.py` | Weekend body authoring — 3 heritage bodies from E-mu XML blueprints. | 2026-04-07 | no | yes — {body.name}.json, trench_live.json |
| `C:/Users/hooki/Trench/trench-juce/forge/pyruntime/recipes/slammed_dark_bright_belch.py` | One usable forge recipe: slammed dark->bright belch. | 2026-04-07 | no | yes — {body.name}.json, {body.name}.report.json |
| `C:/Users/hooki/Trench/trench-juce/tools/author_speaker_knockerz_current_runtime_v1.py` | (no docstring) | 2026-04-07 | no | yes — P2k_006.json, Speaker_Knockerz_current_runtime_v1.keyframe.json |
| `C:/Users/hooki/Trench/trench-juce/tools/author_speaker_knockerz_v3_repro.py` | (no docstring) | 2026-04-07 | no | yes — P2k_006.json, Speaker_Knockerz_v3.keyframe.json |
| `C:/Users/hooki/Trench/trench-juce/tools/forge_pymoo.py` | Multi-objective optimization of TRENCH bodies via pymoo. | 2026-04-07 | no | yes — {name}.json, _summary.json |
| `C:/Users/hooki/Trench/trench-juce/tools/mass_generate.py` | Mass candidate generator for all 4 shipping bodies. | 2026-04-07 | no | yes — _summary.json |
| `C:/Users/hooki/Trench/trench-juce/tools/measure_body.py` | Measure a body against the Speaker Knockerz shipping gates. | 2026-04-07 | no | no |
| `C:/Users/hooki/Trench/trench-juce/tools/plot_body_comparison.py` | Plot frequency responses comparing two bodies at key morph/Q positions. | 2026-04-07 | no | no |

### Rust binaries (2)

| Path | Purpose | Last Modified | Referenced By | Writes Files |
|------|---------|---------------|---------------|--------------|
| `C:/Users/hooki/Trench/trench-juce/trench-mcp/src/main.rs` | ! trench-mcp — MCP server exposing TRENCH analysis to Claude. | 2026-04-11 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/Trench/xtask/src/main.rs` | (no docstring) | 2026-04-12 | yes — Trench\dev\sift\gen_inventory.py | — |

## C:/Users/hooki/trench_re_vault

### Python (33 scripts)

| Path | Purpose | Last Modified | Referenced By | Writes Files |
|------|---------|---------------|---------------|--------------|
| `C:/Users/hooki/trench_re_vault/ghidra_scripts/e4ultra_headless_init.py` | @category E4Ultra | 2026-03-04 | no | no |
| `C:/Users/hooki/trench_re_vault/tools/codec_lab/run_dual_target_dirty.py` | Run G.726 then CVSD DIRTY sessions with fixed ordering. | 2026-03-05 | no | no |
| `C:/Users/hooki/trench_re_vault/tools/cvsd_lab/run_cvsd_end_to_end.py` | DIRTY end-to-end behavioral extraction for CVSD. | 2026-03-05 | yes — trench_re_vault\tools\codec_lab\run_dual_target_dirty.py | yes — CAPTURE_MANIFEST.json, cvsd_behavioral_spec.v1.yaml |
| `C:/Users/hooki/trench_re_vault/tools/e4ultra/extract_morph_designer_tables.py` | extract_morph_designer_tables.py — Extract CPhantomMorphDesigner lookup tables from EmulatorX.dll | 2026-03-08 | no | yes — analysis/morph_designer_tables.json |
| `C:/Users/hooki/trench_re_vault/tools/e4ultra/hash_roms.py` | Hash EOS/E4 Ultra ROM images and emit manifest.json. | 2026-03-04 | no | yes — analysis/e4ultra/roms/manifest.json |
| `C:/Users/hooki/trench_re_vault/tools/g726_lab/run_g726_end_to_end.py` | DIRTY end-to-end behavioral extraction for G.726 ADPCM. | 2026-03-05 | yes — trench_re_vault\tools\codec_lab\run_dual_target_dirty.py | yes — probe_{name}_cs{code_size}.bin, transfer_curve_{mode_tag}.csv |
| `C:/Users/hooki/trench_re_vault/tools/g729_lab/run_g729_decoder_first_extract.py` | DIRTY decoder-first behavioral extraction for a G.729-class LSP grid target. | 2026-03-06 | no | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/analyze_dc_curve.py` | (no docstring) | 2026-03-09 | yes — trench_re_vault\tools\inflator_lab\run_inflator_blackbox_capture.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/analyze_harmonics.py` | (no docstring) | 2026-03-09 | yes — trench_re_vault\tools\inflator_lab\run_inflator_blackbox_capture.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/analyze_memory.py` | (no docstring) | 2026-03-03 | yes — trench_re_vault\tools\inflator_lab\run_inflator_blackbox_capture.py | yes — memory_report.json, memory_report.txt |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/build_inflator_handoff.py` | (no docstring) | 2026-03-07 | no | yes — MANIFEST.json, measurement_manifest.json |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/build_inflator_native_handoff.py` | (no docstring) | 2026-03-09 | no | yes — MANIFEST.json, measurement_manifest.json |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/extract_inflator_re_snapshot.py` | (no docstring) | 2026-03-10 | no | yes |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/generate_all_inputs.py` | (no docstring) | 2026-03-03 | no | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/generate_bandlimited.py` | (no docstring) | 2026-03-03 | yes — trench_re_vault\tools\inflator_lab\generate_all_inputs.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/generate_burst.py` | (no docstring) | 2026-03-03 | yes — trench_re_vault\tools\inflator_lab\generate_all_inputs.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/generate_dc_steps.py` | (no docstring) | 2026-03-03 | yes — trench_re_vault\tools\inflator_lab\generate_all_inputs.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/generate_multitone.py` | (no docstring) | 2026-03-03 | yes — trench_re_vault\tools\inflator_lab\generate_all_inputs.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/generate_sine_levels.py` | (no docstring) | 2026-03-03 | yes — trench_re_vault\tools\inflator_lab\generate_all_inputs.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/generate_swept_sine.py` | (no docstring) | 2026-03-03 | yes — trench_re_vault\tools\inflator_lab\generate_all_inputs.py | no |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/run_inflator_blackbox_capture.py` | (no docstring) | 2026-03-09 | yes — trench_re_vault\tools\inflator_lab\build_inflator_native_handoff.py | yes — sine_manifest_sr44100.json, burst_manifest_sr44100.json |
| `C:/Users/hooki/trench_re_vault/tools/inflator_lab/run_inflator_full_capture.py` | (no docstring) | 2026-03-09 | no | yes — determinism_report.json, host_probe.txt |
| `C:/Users/hooki/trench_re_vault/tools/morpheus_lab/bake_responses.py` | Bake Morpheus cubes into a Trenchwork-ready sandbox dataset. | 2026-03-05 | yes — trenchwork_clean\tools\ai_lab\build_sandbox_session.py | yes — cube_{cube_index:03d}.json, ]:03d}.json |
| `C:/Users/hooki/trench_re_vault/tools/morpheus_lab/decode_cubes.py` | Decode Morpheus cubes from the legacy 8-sample-per-bit WAV container. | 2026-03-04 | no | yes — artifacts/morpheus_cubes_decoded.json |
| `C:/Users/hooki/trench_re_vault/tools/morpheus_lab/plot_p2k_stage_responses_20260312.py` | Plot isolated and cascaded responses for every stage in a P2k source cartridge. | 2026-03-12 | no | yes — C:\Users\hooki\do-it\cartridges\P2k_013.json, CAPTURE_MANIFEST.json |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/build_measure_harness.py` | QSound phase-5 measurement harness scaffolding. | 2026-03-04 | yes — trenchwork_clean\tools\offline-render\src\main.rs | yes — measurement_manifest.json |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/build_phase1_queue.py` | Build a Phase-1 rename/decompile queue from qsound_table_banks JSON output. | 2026-03-03 | no | yes — C:/Users/hooki/Downloads/qsound_table_banks_*.json, tools/qsound_lab/QSOUND_PHASE1_RENAME_QUEUE.md |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/build_qsound_phase3_bundle.py` | Build QSound phase-3 extraction bundle docs from phase1/phase2/decompile artifacts. | 2026-03-03 | no | yes — QSOUND_RENDERER_A_SPEC.md, QSOUND_RENDERER_C_SPEC.md |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/build_qsound_sanitized_bundle.py` | Build a sanitized QSound handoff bundle inside DIRTY and enforce fail-closed leak gates. | 2026-03-07 | yes — trench_re_vault\tools\qsound_lab\run_qsound_spatial_batch_capture.py | yes — HANDOFF.json, CONTRACT.json |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/extract_qcreator.py` | Extract QCreator payload files from legacy InstallShield artifacts. | 2026-03-03 | no | yes — extract_manifest.json |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/reconstruct_qsound_spatial_v1.py` | QSound spatial v1 clean-room reconstruction pipeline. | 2026-03-07 | yes — trench_re_vault\tools\qsound_lab\run_qsound_spatial_batch_capture.py | yes — train_metrics.csv, val_metrics.csv |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/run_qsound_pan_batch_capture.py` | Manifest-driven QSound pan-grid batch capture for DIRTY workspace. | 2026-03-05 | no | yes — Batch-capture QSound pan grid and emit CAPTURE_MANIFEST.json, CAPTURE_PLAN.json |
| `C:/Users/hooki/trench_re_vault/tools/qsound_lab/run_qsound_spatial_batch_capture.py` | QSound live-ingress v2 batch capture driver for the DIRTY airlock. | 2026-03-07 | no | yes — artifacts/qsound_phase5_harness/measurement_manifest.json, stimulus_variants.json |

## C:/Users/hooki/trenchwork_clean

### Python (31 scripts)

| Path | Purpose | Last Modified | Referenced By | Writes Files |
|------|---------|---------------|---------------|--------------|
| `C:/Users/hooki/trenchwork_clean/files/forge_v2.py` | forge_v2.py — Structured body choreography compiler for TRENCH. | 2026-04-03 | no | yes |
| `C:/Users/hooki/trenchwork_clean/pyruntime/forge_optimize.py` | Forge optimizer — pymoo multi-objective body search. | 2026-04-07 | yes — Trench\trench-juce\forge\pyruntime\forge_heritage_optimize.py | yes — {name}.json, pareto.png |
| `C:/Users/hooki/trenchwork_clean/pyruntime/forge_plot.py` | Plot frequency responses for TRENCH bodies. | 2026-04-04 | no | no |
| `C:/Users/hooki/trenchwork_clean/pyruntime/recipes/author_weekend.py` | Weekend body authoring — 3 heritage bodies from E-mu XML blueprints. | 2026-03-28 | no | yes — {body.name}.json, trench_live.json |
| `C:/Users/hooki/trenchwork_clean/pyruntime/recipes/slammed_dark_bright_belch.py` | One usable forge recipe: slammed dark->bright belch. | 2026-03-28 | no | yes — {body.name}.json, {body.name}.report.json |
| `C:/Users/hooki/trenchwork_clean/tools/ai_lab/autonomous_oddness_loop.py` | Autonomous oddness-focused generation loop. | 2026-03-23 | no | yes — manifest.json, manifest.curves.bin |
| `C:/Users/hooki/trenchwork_clean/tools/ai_lab/build_sandbox_session.py` | One-shot sandbox session builder: | 2026-03-23 | no | yes — session_provenance.json, decoded json OR manifest.json |
| `C:/Users/hooki/trenchwork_clean/tools/ai_lab/rank_topk.py` | Rank entries in a manifest and select top-k. | 2026-03-23 | yes — trenchwork_clean\tools\ai_lab\autonomous_oddness_loop.py | yes — topk_manifest.curves.bin, topk_manifest.json |
| `C:/Users/hooki/trenchwork_clean/tools/ai_lab/repair_labels_jsonl.py` | Repair malformed labels JSONL and enforce last-write-wins per id. | 2026-03-23 | no | no |
| `C:/Users/hooki/trenchwork_clean/tools/apply_verdicts.py` | Apply triage verdicts from the browser terminal to the filesystem. | 2026-03-23 | no | yes — site/observatory_data.json, site/observatory/public/observatory_data.json |
| `C:/Users/hooki/trenchwork_clean/tools/calibration_to_body.py` | calibration_to_body.py — Convert calibration JSON directly to compiled-v1 body. | 2026-04-03 | no | yes |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/acceptance_tests.py` | (no docstring) | 2026-03-23 | no | yes |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/build_synthetic_session.py` | First-principles synthetic response family: | 2026-03-23 | no | yes — synthetic/{rid}.json, manifest.curves.bin |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/check_clean_inputs.py` | (no docstring) | 2026-03-23 | yes — trenchwork_clean\tools\cleanroom\run_cleanroom_gate.py | no |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/extract_spatial_handoff.py` | (no docstring) | 2026-03-23 | yes — trenchwork_clean\tools\cleanroom\run_cleanroom_gate.py | yes — CONTRACT.json, CONTRACT.json |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/promote_handoff.py` | (no docstring) | 2026-03-23 | yes — trenchwork_clean\tools\cleanroom\promote_qsound_spatial_v2.py | yes — handoff_bundle.v1.schema.json, handoff_contract.v1.schema.json |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/promote_qsound_spatial_v2.py` | (no docstring) | 2026-03-23 | no | yes — measurement_manifest.json, contracts/cleanroom/qsound_spatial_acceptance.v1.json |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/rank_session.py` | (no docstring) | 2026-03-23 | yes — trenchwork_clean\tools\cleanroom\build_synthetic_session.py | yes — Path to manifest.json, score_overlay.json |
| `C:/Users/hooki/trenchwork_clean/tools/cleanroom/run_cleanroom_gate.py` | (no docstring) | 2026-03-23 | no | no |
| `C:/Users/hooki/trenchwork_clean/tools/emu_historical_forge.py` | emu_historical_forge.py — Compile raw E-mu MorphDesigner templates through | 2026-04-02 | no | yes — Historical_{key}.json |
| `C:/Users/hooki/trenchwork_clean/tools/forge_batch.py` | forge_batch.py — Cancellation Engine body generator for TRENCH Sift workflow. | 2026-04-01 | no | yes — {i:03d}_{name}.json |
| `C:/Users/hooki/trenchwork_clean/tools/forge_generator.py` | forge_generator.py — Anchor + Vector body generator for TRENCH Sift workflow. | 2026-04-03 | yes — Trench\trench-juce\forge\pyruntime\forge_shipping.py | yes |
| `C:/Users/hooki/trenchwork_clean/tools/measure_x3_parity.py` | Measure spectral parity between X3 reference captures and TRENCH offline renders. | 2026-03-23 | no | yes |
| `C:/Users/hooki/trenchwork_clean/tools/metric_audit_path.py` | Summarize benchmark-vs-reference audit support from repo assets. | 2026-03-23 | no | yes — p2k_filter_names.json, CONTRACT.json |
| `C:/Users/hooki/trenchwork_clean/tools/oasys_lab/extract_oasys_payloads.py` | Deterministic extraction + triage for legacy OASYS-related archive drops. | 2026-03-23 | no | yes — oasys_extract_report_{ts}.json, oasys_extract_report_{ts}.md |
| `C:/Users/hooki/trenchwork_clean/tools/plot_parity.py` | Plot TRENCH vs X3 transfer function comparison at all 4 corners. | 2026-03-23 | no | no |
| `C:/Users/hooki/trenchwork_clean/tools/qsound_spatialize.py` | QSound HRIR spatializer — post-body spatial placement tool. | 2026-03-29 | no | no |
| `C:/Users/hooki/trenchwork_clean/tools/render_crawl.py` | Render Space Crawler CSV grids as heatmap PNGs. | 2026-03-23 | no | no |
| `C:/Users/hooki/trenchwork_clean/tools/render_sample_rate_comparison.py` | Render a raw-capture cartridge response at two sample rates. | 2026-03-23 | no | yes — {args.slug}.png, {args.slug}.json |
| `C:/Users/hooki/trenchwork_clean/tools/sift.py` | sift.py — Tinder for DSP. | 2026-04-01 | yes — Trench\dev\sift\gen_inventory.py | no |
| `C:/Users/hooki/trenchwork_clean/triage.py` | Repo-root shim to the active sift workflow. | 2026-04-04 | yes — trenchwork_clean\tools\apply_verdicts.py | no |

### Rust binaries (10)

| Path | Purpose | Last Modified | Referenced By | Writes Files |
|------|---------|---------------|---------------|--------------|
| `C:/Users/hooki/trenchwork_clean/archive/pre-runtime/trench-station/src/main.rs` | ! TRENCH Station — standalone workstation for body audition. | 2026-03-26 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/archive/tools/migrate-atlas/src/main.rs` | (no docstring) | 2026-03-23 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/foo/src/main.rs` | (no docstring) | 2026-03-27 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/quarantine/trench-ui/src/main.rs` | Resolve workspace root — works whether CWD is workspace root or trench-ui/. | 2026-03-23 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/tools/bake/src/main.rs` | ! One-time bake tool: legacy cartridge JSON → compiled JSON + const Rust source. | 2026-03-23 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/tools/offline-render/src/main.rs` | Render QSound stimuli through a simplified clean-room path using extracted tables. | 2026-03-23 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/tools/render-morph/src/main.rs` | QSound contract ITD/ILD coefficients (from CONTRACT.json behavioral_laws_v1) | 2026-03-23 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/tools/render-parity/src/main.rs` | ! Render P2K_013 (Talking Hedz) at all 4 corners through the TRENCH cascade. | 2026-03-23 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/trench-mcp/src/main.rs` | ! trench-mcp — MCP server exposing TRENCH analysis to Claude. | 2026-03-25 | yes — Trench\dev\sift\gen_inventory.py | — |
| `C:/Users/hooki/trenchwork_clean/trench-plugin/src/main.rs` | WASAPI may not honor the default 512 period size and can deliver more | 2026-04-04 | yes — Trench\dev\sift\gen_inventory.py | — |

## Duplicate stems

- **author_weekend** (2x): C:/Users/hooki/Trench/trench-juce/forge/pyruntime/recipes/author_weekend.py | C:/Users/hooki/trenchwork_clean/pyruntime/recipes/author_weekend.py
- **calibration_to_body** (2x): C:/Users/hooki/Trench/trench-juce/forge/pyruntime/calibration_to_body.py | C:/Users/hooki/trenchwork_clean/tools/calibration_to_body.py
- **forge_generator** (2x): C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_generator.py | C:/Users/hooki/trenchwork_clean/tools/forge_generator.py
- **forge_optimize** (2x): C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_optimize.py | C:/Users/hooki/trenchwork_clean/pyruntime/forge_optimize.py
- **forge_plot** (2x): C:/Users/hooki/Trench/trench-juce/forge/pyruntime/forge_plot.py | C:/Users/hooki/trenchwork_clean/pyruntime/forge_plot.py
- **slammed_dark_bright_belch** (2x): C:/Users/hooki/Trench/trench-juce/forge/pyruntime/recipes/slammed_dark_bright_belch.py | C:/Users/hooki/trenchwork_clean/pyruntime/recipes/slammed_dark_bright_belch.py