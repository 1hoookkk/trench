//! Smoke gate: every cartridge under `cartridges/engine/` must load
//! through `Cartridge::from_json` without error.
//!
//! This is intentionally the dumbest useful test — it does not try to
//! measure audio, enforce invariants, or judge taste. It only proves
//! that the shipping phoneme pills baked by `tools/bake_phoneme_pills.py`
//! are all structurally valid `compiled-v1` cartridges the runtime can
//! consume. If this gate fails, either the bake is emitting bad files
//! or the upstream shape bank has drifted.

use std::fs;
use std::path::{Path, PathBuf};

use trench_core::Cartridge;

fn engine_root() -> PathBuf {
    // Tests run with CWD = trench-core/, so cartridges/engine/ is up one level.
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace parent")
        .join("cartridges")
        .join("engine")
}

fn collect_cartridges(dir: &Path, out: &mut Vec<PathBuf>) {
    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return,
    };
    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_cartridges(&path, out);
        } else if path.extension().and_then(|s| s.to_str()) == Some("json")
            && path.file_name().and_then(|s| s.to_str()) != Some("manifest.json")
        {
            out.push(path);
        }
    }
}

#[test]
fn every_engine_pill_loads() {
    let root = engine_root();
    if !root.exists() {
        // Engine dir hasn't been baked yet — treat as skip rather than fail.
        // Run `python tools/bake_phoneme_pills.py` from the repo root to
        // populate it, then re-run this test.
        eprintln!("skip: {} does not exist", root.display());
        return;
    }

    let mut pills = Vec::new();
    collect_cartridges(&root, &mut pills);
    assert!(
        !pills.is_empty(),
        "engine dir {} exists but has no cartridges",
        root.display()
    );

    let mut failures: Vec<String> = Vec::new();
    for pill in &pills {
        let json = match fs::read_to_string(pill) {
            Ok(j) => j,
            Err(e) => {
                failures.push(format!("{}: read error: {e}", pill.display()));
                continue;
            }
        };
        if let Err(e) = Cartridge::from_json(&json) {
            failures.push(format!("{}: {e}", pill.display()));
        }
    }

    assert!(
        failures.is_empty(),
        "{} / {} pills failed to load:\n  {}",
        failures.len(),
        pills.len(),
        failures.join("\n  ")
    );

    println!("loaded {} phoneme pills from {}", pills.len(), root.display());
}
