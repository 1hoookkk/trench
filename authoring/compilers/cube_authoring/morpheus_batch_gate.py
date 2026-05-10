"""Batch-run edge_gate over the morpheus _morpheus/ staging cubes.

Triage scope: `hf_resonance` (sine_sweep, reese, drum_loop, noise) — the
only doctrine-recognized scope whose probes are all present pre-vocal.wav.
This is a pre-screen, not a broadband admission. Cross-references the
trench_re_vault `morpheus_cubes_baked/manifest.json` features so the user
can rank survivors by anomaly_score / centroid_span / dynamic range.

Outputs:
  - cartridges/cubes/_morpheus/_verdicts/<cube>.verdict.json   (one per cube)
  - cartridges/cubes/_morpheus/_verdicts/summary.csv           (one row per cube)

Re-runnable: skips cubes whose verdict file already exists unless --force.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parents[2]
STAGING_DIR = REPO_ROOT / "cartridges" / "cubes" / "_morpheus"
VERDICTS_DIR = STAGING_DIR / "_verdicts"
MANIFEST_PATH = REPO_ROOT.parent / "trench_re_vault" / "artifacts" / "morpheus_cubes_baked" / "manifest.json"

sys.path.insert(0, str(REPO_ROOT / "tools" / "cube_authoring"))
from edge_gate import run_edge_gate  # noqa: E402

TRIAGE_SCOPE = "hf_resonance"


def load_manifest_features(path: Path) -> dict[str, dict]:
    if not path.exists():
        print(f"warn: manifest not found at {path}; running without feature cross-ref", file=sys.stderr)
        return {}
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    return {e["id"]: e for e in data.get("entries", [])}


def summarize_verdict(verdict: dict) -> dict:
    per_probe = verdict.get("per_probe", [])
    fail_counts: dict[str, int] = {}
    edges_passed_total = 0
    edges_total = 0
    for p in per_probe:
        if p.get("status") != "passed" and p.get("status") != "failed":
            continue
        edges_passed_total += p.get("edges_passed", 0)
        edges_total += p.get("edge_count", 0)
        for edge in p.get("edges", []):
            for f in edge.get("fails", []):
                fail_counts[f] = fail_counts.get(f, 0) + 1
    return {
        "edges_passed_total": edges_passed_total,
        "edges_total": edges_total,
        "fail_counts": fail_counts,
        "missing": verdict.get("missing_samples", []),
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--force", action="store_true", help="overwrite existing verdicts")
    parser.add_argument("--scope", default=TRIAGE_SCOPE)
    args = parser.parse_args(argv)

    VERDICTS_DIR.mkdir(parents=True, exist_ok=True)
    cube_paths = sorted(STAGING_DIR.glob("morph_cube_*.cube.json"))
    if not cube_paths:
        print(f"no cubes found in {STAGING_DIR}", file=sys.stderr)
        return 1
    end = len(cube_paths) if args.limit is None else min(len(cube_paths), args.start + args.limit)
    cube_paths = cube_paths[args.start:end]
    features = load_manifest_features(MANIFEST_PATH)

    csv_path = VERDICTS_DIR / "summary.csv"
    csv_exists = csv_path.exists()
    fieldnames = [
        "cube_id", "admitted", "edges_passed", "edges_total",
        "fail_dead", "fail_blow", "fail_collapse", "fail_samey",
        "anomaly_score", "centroid_span_hz", "peak_dynamic_range_db",
        "rms_dynamic_range_db", "terrain_entropy", "peak_hotspot_count",
        "elapsed_s",
    ]
    csv_mode = "a" if csv_exists and not args.force else "w"
    with csv_path.open(csv_mode, newline="", encoding="utf-8") as csvfh:
        writer = csv.DictWriter(csvfh, fieldnames=fieldnames)
        if csv_mode == "w":
            writer.writeheader()

        for i, cp in enumerate(cube_paths):
            cube_id = cp.stem.replace(".cube", "")
            verdict_path = VERDICTS_DIR / f"{cube_id}.verdict.json"
            if verdict_path.exists() and not args.force:
                print(f"[{i+1}/{len(cube_paths)}] skip {cube_id} (verdict exists)")
                continue
            with cp.open(encoding="utf-8") as fh:
                cube = json.load(fh)
            t0 = time.time()
            verdict = run_edge_gate(cube, scope_override=args.scope)
            elapsed = time.time() - t0
            verdict_path.write_text(json.dumps(verdict, indent=2) + "\n", encoding="utf-8")

            summary = summarize_verdict(verdict)
            feat = features.get(cube_id, {})
            row = {
                "cube_id": cube_id,
                "admitted": verdict.get("admitted"),
                "edges_passed": summary["edges_passed_total"],
                "edges_total": summary["edges_total"],
                "fail_dead": summary["fail_counts"].get("dead_zone", 0),
                "fail_blow": summary["fail_counts"].get("blow_up", 0),
                "fail_collapse": summary["fail_counts"].get("abrupt_collapse", 0),
                "fail_samey": summary["fail_counts"].get("samey", 0),
                "anomaly_score": feat.get("anomaly_score"),
                "centroid_span_hz": feat.get("centroid_span_hz"),
                "peak_dynamic_range_db": feat.get("peak_dynamic_range_db"),
                "rms_dynamic_range_db": feat.get("rms_dynamic_range_db"),
                "terrain_entropy": feat.get("terrain_entropy"),
                "peak_hotspot_count": feat.get("peak_hotspot_count"),
                "elapsed_s": round(elapsed, 2),
            }
            writer.writerow(row)
            csvfh.flush()
            print(f"[{i+1}/{len(cube_paths)}] {cube_id} admitted={verdict.get('admitted')} "
                   f"edges={summary['edges_passed_total']}/{summary['edges_total']} "
                   f"fails={summary['fail_counts']} {elapsed:.1f}s")
    print(f"done. csv → {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
