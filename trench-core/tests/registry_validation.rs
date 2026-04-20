//! Validates `FILTER_ARTIFACTS.json` — the single source of truth for
//! filter-artifact identity. The registry is the counter-measure to the
//! display-name collision bug that shipped a false-premise parity test
//! (`hedz_rom` vs raw `P2k_013`) as if they were the same filter.
//!
//! Every record must:
//!   - have a unique, code-safe `internal_id`
//!   - name a concrete repo path that exists on disk
//!   - declare a `derived_from` that resolves to another registered ID
//!     (or is `null` for root sources)
//!   - use only the known enum values for `source_class` and
//!     `truth_status`
//!
//! Every `null_test_pair` and `forbidden_comparison` must name registered
//! IDs. This is how the registry stops being inert documentation.

use std::collections::HashSet;
use std::path::PathBuf;

use serde_json::Value;

fn repo_root() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .expect("workspace parent")
        .to_path_buf()
}

fn load_registry() -> Value {
    let path = repo_root().join("FILTER_ARTIFACTS.json");
    let s = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("read {}: {e}", path.display()));
    serde_json::from_str(&s)
        .unwrap_or_else(|e| panic!("parse FILTER_ARTIFACTS.json: {e}"))
}

/// Strip a `#Lnnn` suffix when resolving a registry path to disk.
fn strip_line_anchor(p: &str) -> &str {
    p.split('#').next().unwrap_or(p)
}

/// `{M0_Q0,…}` glob expansion — the registry uses this for WAV sets.
fn expand_brace(p: &str) -> Vec<String> {
    let (pre, rest) = match p.split_once('{') {
        Some(pair) => pair,
        None => return vec![p.to_string()],
    };
    let (body, post) = match rest.split_once('}') {
        Some(pair) => pair,
        None => return vec![p.to_string()],
    };
    body.split(',')
        .flat_map(|part| expand_brace(&format!("{pre}{part}{post}")))
        .collect()
}

const KNOWN_SOURCE_CLASSES: &[&str] = &[
    "heritage_designer_xml",
    "calibration_re",
    "raw_p2k_skin",
    "compiled_v1",
    "rust_const_rom",
    "canonical_render",
];

const KNOWN_TRUTH_STATUSES: &[&str] = &[
    "authoring_source",
    "runtime_truth",
    "runtime_compile",
    "theoretical_reference",
    "reference_render",
    "deprecated",
];

const KNOWN_Q_BEHAVIORS: &[&str] = &["morph_only", "live_q", "unknown"];

#[test]
fn registry_parses_and_has_schema() {
    let reg = load_registry();
    assert!(reg.get("_schema").is_some(), "missing _schema block");
    assert!(reg["artifacts"].is_array(), "artifacts must be array");
    assert!(reg["null_test_pairs"].is_array(), "null_test_pairs must be array");
    assert!(reg["forbidden_comparisons"].is_array(), "forbidden_comparisons must be array");
}

#[test]
fn registry_ids_are_unique() {
    let reg = load_registry();
    let mut seen = HashSet::new();
    for a in reg["artifacts"].as_array().unwrap() {
        let id = a["internal_id"].as_str().expect("internal_id is string");
        assert!(seen.insert(id.to_string()), "duplicate internal_id: {id}");
    }
}

#[test]
fn registry_enums_are_valid() {
    let reg = load_registry();
    for a in reg["artifacts"].as_array().unwrap() {
        let id = a["internal_id"].as_str().unwrap();
        let sc = a["source_class"].as_str().expect("source_class is string");
        let ts = a["truth_status"].as_str().expect("truth_status is string");
        let qb = a["q_behavior"].as_str().expect("q_behavior is string");
        assert!(
            KNOWN_SOURCE_CLASSES.contains(&sc),
            "{id}: unknown source_class {sc:?}"
        );
        assert!(
            KNOWN_TRUTH_STATUSES.contains(&ts),
            "{id}: unknown truth_status {ts:?}"
        );
        assert!(
            KNOWN_Q_BEHAVIORS.contains(&qb),
            "{id}: unknown q_behavior {qb:?}"
        );
    }
}

#[test]
fn registry_paths_exist_on_disk() {
    let reg = load_registry();
    let root = repo_root();
    for a in reg["artifacts"].as_array().unwrap() {
        let id = a["internal_id"].as_str().unwrap();
        let path_field = a["path"].as_str().expect("path is string");
        let cleaned = strip_line_anchor(path_field);
        for concrete in expand_brace(cleaned) {
            let abs = root.join(&concrete);
            assert!(
                abs.exists(),
                "{id}: path does not exist: {}",
                abs.display()
            );
        }
    }
}

#[test]
fn registry_derived_from_resolves() {
    let reg = load_registry();
    let ids: HashSet<String> = reg["artifacts"]
        .as_array()
        .unwrap()
        .iter()
        .map(|a| a["internal_id"].as_str().unwrap().to_string())
        .collect();
    for a in reg["artifacts"].as_array().unwrap() {
        let id = a["internal_id"].as_str().unwrap();
        match &a["derived_from"] {
            Value::Null => {}
            Value::String(s) => assert!(
                ids.contains(s),
                "{id}: derived_from {s:?} is not a registered internal_id"
            ),
            other => panic!("{id}: derived_from must be string or null, got {other:?}"),
        }
    }
}

#[test]
fn null_test_pairs_reference_registered_ids() {
    let reg = load_registry();
    let ids: HashSet<String> = reg["artifacts"]
        .as_array()
        .unwrap()
        .iter()
        .map(|a| a["internal_id"].as_str().unwrap().to_string())
        .collect();
    for pair in reg["null_test_pairs"].as_array().unwrap() {
        for key in ["predicted", "reference"] {
            let v = pair[key].as_str().expect("pair field is string");
            assert!(
                ids.contains(v),
                "null_test_pair.{key} = {v:?} not in registry"
            );
        }
    }
    for pair in reg["forbidden_comparisons"].as_array().unwrap() {
        for key in ["predicted", "reference"] {
            let v = pair[key].as_str().expect("pair field is string");
            assert!(
                ids.contains(v),
                "forbidden_comparison.{key} = {v:?} not in registry"
            );
        }
    }
}
