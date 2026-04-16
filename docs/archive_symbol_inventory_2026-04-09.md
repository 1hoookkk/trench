# Archive Symbol Inventory

- Archive root: `C:/Users/hooki/trenchwork_clean/archive`
- Live root: `C:/Users/hooki/Trench`
- Total symbols: `1134`
- `alive_in_shipping`: `50`
- `dead_heritage`: `98`
- `portable_to_current`: `863`
- `needs_adaptation`: `146`

## Focus Symbol Status

| Symbol | Status | Archive Source | Live Hits |
|---|---|---|---|
| InterferenceLink | portable_to_current | pre-runtime/trench-forge/src/cluster.rs:56 | - |
| ClusterPlan | portable_to_current | pre-runtime/trench-forge/src/cluster.rs:69 | - |
| Role | portable_to_current | pre-runtime/trench-forge/src/coordination.rs:38 | - |
| owned_gaps | portable_to_current | pre-runtime/trench-forge/src/coordination.rs:399 | - |
| seed_relative_zero | portable_to_current | pre-runtime/trench-forge/src/coordination.rs:481 | - |
| pressure_law | portable_to_current | pre-runtime/trench-forge/src/redistribution.rs:145 | - |
| prominence | portable_to_current | pre-runtime/trench-forge/src/redistribution.rs:161 | - |
| topology_anchor_theft | portable_to_current | pre-runtime/trench-forge/src/coordination.rs:1548 | - |
| topology_interference_drift | portable_to_current | pre-runtime/trench-forge/src/coordination.rs:1587 | - |
| topology_cluster_leakage | portable_to_current | pre-runtime/trench-forge/src/coordination.rs:1626 | - |
| topology_bypass_wake | portable_to_current | pre-runtime/trench-forge/src/coordination.rs:1666 | - |
| Material | needs_adaptation | trench-atlas/src/material.rs:140 | - |
| Motor | needs_adaptation | trench-atlas/src/motor.rs:17 | - |
| Volatility | needs_adaptation | trench-atlas/src/material.rs:75 | - |
| Curation | needs_adaptation | trench-atlas/src/material.rs:63 | - |
| House | needs_adaptation | trench-atlas/src/material.rs:43 | - |

## Portable Candidates (Top 50)

| Symbol | Kind | Source |
|---|---|---|
| StageFunction | pub_enum | pre-runtime/trench-forge/src/cluster.rs:12 |
| ClusterStage | pub_struct | pre-runtime/trench-forge/src/cluster.rs:25 |
| Cluster | pub_struct | pre-runtime/trench-forge/src/cluster.rs:33 |
| CompiledTopology | pub_struct | pre-runtime/trench-forge/src/cluster.rs:41 |
| InterferenceLink | pub_struct | pre-runtime/trench-forge/src/cluster.rs:56 |
| ClusterPlan | pub_struct | pre-runtime/trench-forge/src/cluster.rs:69 |
| active_count | pub_fn | pre-runtime/trench-forge/src/cluster.rs:77 |
| subordinate_count | pub_fn | pre-runtime/trench-forge/src/cluster.rs:85 |
| total_stages | pub_fn | pre-runtime/trench-forge/src/cluster.rs:93 |
| resolve_links | fn | pre-runtime/trench-forge/src/cluster.rs:167 |
| to_skeleton | pub_fn | pre-runtime/trench-forge/src/cluster.rs:179 |
| seed_zeros | pub_fn | pre-runtime/trench-forge/src/cluster.rs:184 |
| FourCornerTopology | pub_struct | pre-runtime/trench-forge/src/cluster.rs:318 |
| build_four_corners | pub_fn | pre-runtime/trench-forge/src/cluster.rs:331 |
| inherit_q | fn | pre-runtime/trench-forge/src/cluster.rs:351 |
| check_q_inheritance | pub_fn | pre-runtime/trench-forge/src/cluster.rs:404 |
| check_pair | fn | pre-runtime/trench-forge/src/cluster.rs:412 |
| default_cluster_layout | pub_fn | pre-runtime/trench-forge/src/cluster.rs:439 |
| four_pole_plan | fn | pre-runtime/trench-forge/src/cluster.rs:474 |
| six_pole_plan | fn | pre-runtime/trench-forge/src/cluster.rs:483 |
| default_layout_valid | fn | pre-runtime/trench-forge/src/cluster.rs:495 |
| default_layout_six_poles | fn | pre-runtime/trench-forge/src/cluster.rs:503 |
| to_skeleton_produces_correct_count | fn | pre-runtime/trench-forge/src/cluster.rs:511 |
| to_skeleton_six_poles | fn | pre-runtime/trench-forge/src/cluster.rs:520 |
| seed_zeros_produces_seeds | fn | pre-runtime/trench-forge/src/cluster.rs:527 |
| bypass_stages_excluded_from_skeleton | fn | pre-runtime/trench-forge/src/cluster.rs:539 |
| too_many_stages_flagged | fn | pre-runtime/trench-forge/src/cluster.rs:563 |
| interference_pair_seeding | fn | pre-runtime/trench-forge/src/cluster.rs:572 |
| unlinked_interference_flagged | fn | pre-runtime/trench-forge/src/cluster.rs:622 |
| self_link_flagged | fn | pre-runtime/trench-forge/src/cluster.rs:642 |
| compile_returns_all_three | fn | pre-runtime/trench-forge/src/cluster.rs:663 |
| compile_preserves_function_tags | fn | pre-runtime/trench-forge/src/cluster.rs:684 |
| q_corners_inherit_topology | fn | pre-runtime/trench-forge/src/cluster.rs:713 |
| q_corners_preserve_freq_vary_radius | fn | pre-runtime/trench-forge/src/cluster.rs:728 |
| q_anchor_radius_preserved | fn | pre-runtime/trench-forge/src/cluster.rs:753 |
| full_pipeline_cluster_to_solver | fn | pre-runtime/trench-forge/src/cluster.rs:769 |
| Role | pub_enum | pre-runtime/trench-forge/src/coordination.rs:38 |
| classify | fn | pre-runtime/trench-forge/src/coordination.rs:49 |
| ZeroStats | pub_struct | pre-runtime/trench-forge/src/coordination.rs:65 |
| CorpusPatterns | pub_struct | pre-runtime/trench-forge/src/coordination.rs:76 |
| analyze_corpus | pub_fn | pre-runtime/trench-forge/src/coordination.rs:85 |
| MidpointReport | pub_struct | pre-runtime/trench-forge/src/coordination.rs:157 |
| probe_midpoints | pub_fn | pre-runtime/trench-forge/src/coordination.rs:174 |
| MidpointAudit | pub_struct | pre-runtime/trench-forge/src/coordination.rs:216 |
| midpoint_gate | pub_fn | pre-runtime/trench-forge/src/coordination.rs:242 |
| SkeletonPlan | pub_struct | pre-runtime/trench-forge/src/coordination.rs:341 |
| CoordinatedCorner | pub_struct | pre-runtime/trench-forge/src/coordination.rs:348 |
| RelativeZero | pub_struct | pre-runtime/trench-forge/src/coordination.rs:361 |
| to_absolute | pub_fn | pre-runtime/trench-forge/src/coordination.rs:367 |
| OwnedGap | pub_struct | pre-runtime/trench-forge/src/coordination.rs:380 |

## Needs Adaptation (Top 50)

| Symbol | Kind | Source |
|---|---|---|
| SawOsc | type_item | pre-runtime/trench-station/src/audio.rs:16 |
| set_freq | fn | pre-runtime/trench-station/src/audio.rs:31 |
| next | fn | pre-runtime/trench-station/src/audio.rs:35 |
| SuperSawOsc | type_item | pre-runtime/trench-station/src/audio.rs:46 |
| next | fn | pre-runtime/trench-station/src/audio.rs:63 |
| WavCursor | type_item | pre-runtime/trench-station/src/audio.rs:86 |
| trigger | fn | pre-runtime/trench-station/src/audio.rs:101 |
| next | fn | pre-runtime/trench-station/src/audio.rs:107 |
| RackPlayback | type_item | pre-runtime/trench-station/src/audio.rs:126 |
| refresh_buffers | fn | pre-runtime/trench-station/src/audio.rs:151 |
| continuous_sample | fn | pre-runtime/trench-station/src/audio.rs:161 |
| oneshot_sample | fn | pre-runtime/trench-station/src/audio.rs:193 |
| print_audio_info | pub_fn | pre-runtime/trench-station/src/audio.rs:217 |
| start_audio | pub_fn | pre-runtime/trench-station/src/audio.rs:242 |
| process_audio | fn | pre-runtime/trench-station/src/audio.rs:358 |
| BodyEntry | pub_struct | pre-runtime/trench-station/src/body_loader.rs:15 |
| BodySource | pub_enum | pre-runtime/trench-station/src/body_loader.rs:22 |
| list_bodies | pub_fn | pre-runtime/trench-station/src/body_loader.rs:29 |
| cartridge_to_corners | fn | pre-runtime/trench-station/src/body_loader.rs:120 |
| corners_to_baked | pub_fn | pre-runtime/trench-station/src/body_loader.rs:136 |
| run | pub_fn | pre-runtime/trench-station/src/cli.rs:10 |
| usage | fn | pre-runtime/trench-station/src/cli.rs:23 |
| cmd_list | fn | pre-runtime/trench-station/src/cli.rs:38 |
| cmd_sweep_all | fn | pre-runtime/trench-station/src/cli.rs:126 |
| Xor64 | type_item | pre-runtime/trench-station/src/panels/author.rs:31 |
| next | fn | pre-runtime/trench-station/src/panels/author.rs:34 |
| uniform | fn | pre-runtime/trench-station/src/panels/author.rs:42 |
| roll_author | fn | pre-runtime/trench-station/src/panels/author.rs:51 |
| AuthoredStage | pub_struct | pre-runtime/trench-station/src/panels/author.rs:119 |
| CaptureState | pub_struct | pre-runtime/trench-station/src/panels/author.rs:138 |
| boost_radius | fn | pre-runtime/trench-station/src/panels/author.rs:166 |
| build_from_captures | pub_fn | pre-runtime/trench-station/src/panels/author.rs:187 |
| AuthorState | pub_struct | pre-runtime/trench-station/src/panels/author.rs:205 |
| fork_from | pub_fn | pre-runtime/trench-station/src/panels/author.rs:257 |
| show | pub_fn | pre-runtime/trench-station/src/panels/author.rs:340 |
| show_capture | pub_fn | pre-runtime/trench-station/src/panels/author.rs:419 |
| show_save | pub_fn | pre-runtime/trench-station/src/panels/author.rs:499 |
| BrowserState | pub_struct | pre-runtime/trench-station/src/panels/browser.rs:7 |
| reload | pub_fn | pre-runtime/trench-station/src/panels/browser.rs:20 |
| PromoteResult | pub_struct | pre-runtime/trench-station/src/panels/browser.rs:31 |
| show | pub_fn | pre-runtime/trench-station/src/panels/browser.rs:37 |
| show | pub_fn | pre-runtime/trench-station/src/panels/canvas.rs:12 |
| draw_grid | fn | pre-runtime/trench-station/src/panels/canvas.rs:38 |
| draw_response | fn | pre-runtime/trench-station/src/panels/canvas.rs:81 |
| freq_to_x | fn | pre-runtime/trench-station/src/panels/canvas.rs:108 |
| db_to_y | fn | pre-runtime/trench-station/src/panels/canvas.rs:113 |
| show | pub_fn | pre-runtime/trench-station/src/panels/drumpad.rs:18 |
| author | pub_item | pre-runtime/trench-station/src/panels/mod.rs:1 |
| browser | pub_item | pre-runtime/trench-station/src/panels/mod.rs:2 |
| canvas | pub_item | pre-runtime/trench-station/src/panels/mod.rs:3 |

