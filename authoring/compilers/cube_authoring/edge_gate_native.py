"""Pole-native edge gate for morpheus cubes.

Mirror of `edge_gate.py`, but interpolates **(freq_hz, radius)** along each
edge instead of biquad coefficients. At each of the 5 edge positions the
poles are linearly blended, then converted to a fresh biquad cascade with
passthrough numerator (b0=1, b1=0, b2=0). This is how the Morpheus's
hardware actually swept its cube — it tells us whether the *layouts* are
musical, separate from whether our coeff-interpolating runtime can carry
them between corners.

Verdict field is named `native_edge_pass` and the schema is
`trench.gate.cube_verdict_native.v1` so this is never confused with the
runtime (coeff-space) verdict.

Source of truth: `trench_re_vault/artifacts/morpheus_cubes_decoded.json`.
No compiled-cube staging dependency — operates directly on pole bytes.
"""
from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from itertools import combinations
from pathlib import Path

import numpy as np

REPO_ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO_ROOT / "tools" / "cube_authoring"))
sys.path.insert(0, str(REPO_ROOT / "tools"))

from pole_math import pole_to_kernel_stage  # noqa: E402
from probes import (MissingProbeSample, available_probes, get_probe,  # noqa: E402
                    in_scope_probes)
from edge_gate import (ABRUPT_COLLAPSE_FRACTION, BLOW_UP_DB, DEAD_ZONE_DB,  # noqa: E402
                        EDGE_POSITIONS, EDGES, SAMEY_CENTROID_OCTAVE,
                        SAMEY_L2_MULT, l2_spectrum, level_normalize,
                        magnitude_db, mean_pair_l2, octave_shift,
                        render_probe, spectral_centroid)

CORNER_LABELS = ("c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111")
RUNTIME_SR = 44100
TRIAGE_SCOPE = "hf_resonance"
DEFAULT_VAULT = REPO_ROOT.parent / "trench_re_vault" / "artifacts" / "morpheus_cubes_decoded.json"
DEFAULT_OUT = REPO_ROOT / "cartridges" / "cubes" / "_morpheus" / "_verdicts_native"
DEFAULT_MANIFEST = REPO_ROOT.parent / "trench_re_vault" / "artifacts" / "morpheus_cubes_baked" / "manifest.json"


# -----------------------------------------------------------------------------
# Pole → coeffs (passthrough numerator)
# -----------------------------------------------------------------------------

def pole_to_stage(freq_hz: float, radius: float, sr: int
                   ) -> tuple[float, float, float, float, float]:
    return pole_to_kernel_stage(freq_hz, radius, sr)


def corner_pole_pairs(corner_in: list[dict]) -> list[tuple[float, float]]:
    return [(float(s["freq_hz"]), float(s["radius"])) for s in corner_in]


def interp_poles(a: list[tuple[float, float]], b: list[tuple[float, float]],
                  t: float) -> list[tuple[float, float]]:
    return [(av[0] * (1 - t) + bv[0] * t, av[1] * (1 - t) + bv[1] * t)
             for av, bv in zip(a, b)]


def poles_to_stages(poles: list[tuple[float, float]], sr: int
                     ) -> list[tuple[float, float, float, float, float]]:
    return [pole_to_stage(f, r, sr) for f, r in poles]


# -----------------------------------------------------------------------------
# Per-edge gate (pole-domain interpolation)
# -----------------------------------------------------------------------------

def render_corner_pole(poles: list[tuple[float, float]], probe: np.ndarray,
                        sr: int) -> np.ndarray:
    stages = poles_to_stages(poles, sr)
    out = render_probe(stages, 1.0, probe)
    return level_normalize(out, probe)


def evaluate_edge_native(corner_a_poles, corner_b_poles, probe, probe_max_db,
                          mean_pair, sr) -> dict:
    freqs = None
    spectra_db = []
    max_bins = []
    for t in EDGE_POSITIONS:
        poles_t = interp_poles(corner_a_poles, corner_b_poles, t)
        stages = poles_to_stages(poles_t, sr)
        out = render_probe(stages, 1.0, probe)
        out = level_normalize(out, probe)
        fr, db = magnitude_db(out, sr)
        freqs = fr
        spectra_db.append(db)
        max_bins.append(float(db.max()))
    spectra_db = np.array(spectra_db)

    dead_positions = [v for v in max_bins if v < probe_max_db + DEAD_ZONE_DB]
    blown_positions = [v for v in max_bins if v > probe_max_db + BLOW_UP_DB]
    endpoint_l2 = l2_spectrum(spectra_db[0], spectra_db[-1])
    centroid_a = spectral_centroid(freqs, spectra_db[0])
    centroid_b = spectral_centroid(freqs, spectra_db[-1])
    centroid_oct = octave_shift(centroid_a, centroid_b)
    step_l2s = [l2_spectrum(spectra_db[i], spectra_db[i + 1]) for i in range(4)]
    max_step = max(step_l2s) if step_l2s else 0.0
    collapse_frac = (max_step / endpoint_l2) if endpoint_l2 > 1e-9 else 0.0
    samey = (endpoint_l2 < SAMEY_L2_MULT * mean_pair) and (centroid_oct < SAMEY_CENTROID_OCTAVE)

    fails = []
    if dead_positions:
        fails.append("dead_zone")
    if blown_positions:
        fails.append("blow_up")
    if collapse_frac > ABRUPT_COLLAPSE_FRACTION:
        fails.append("abrupt_collapse")
    if samey:
        fails.append("samey")
    return {"passed": not fails, "fails": fails,
            "endpoint_l2": round(endpoint_l2, 4),
            "centroid_octaves": round(centroid_oct, 4),
            "abrupt_collapse_fraction": round(collapse_frac, 4),
            "max_bin_db_per_position": [round(v, 2) for v in max_bins]}


def run_probe_gate_native(corners_poles: dict[str, list[tuple[float, float]]],
                           probe: np.ndarray, probe_id: str, sr: int) -> dict:
    probe_max_db = float(magnitude_db(probe, sr)[1].max())
    spectra = {label: magnitude_db(render_corner_pole(p, probe, sr), sr)[1]
                for label, p in corners_poles.items()}
    mean_pair = mean_pair_l2(spectra)
    edge_results = []
    for a, b in EDGES:
        r = evaluate_edge_native(corners_poles[a], corners_poles[b], probe,
                                   probe_max_db, mean_pair, sr)
        r["edge"] = f"{a}-{b}"
        edge_results.append(r)
    all_pass = all(r["passed"] for r in edge_results)
    return {"probe": probe_id, "passed": all_pass,
            "edges_passed": sum(1 for r in edge_results if r["passed"]),
            "edge_count": len(edge_results),
            "probe_max_bin_db": round(probe_max_db, 2),
            "mean_cube_pair_l2": round(mean_pair, 4),
            "edges": edge_results}


def run_native_edge_gate(cube_decoded: dict, *, sr: int = RUNTIME_SR,
                           scope: str = TRIAGE_SCOPE) -> dict:
    corners_in = cube_decoded["corners"]
    corners_poles = {label: corner_pole_pairs(corners_in[i])
                       for i, label in enumerate(CORNER_LABELS)}
    probe_ids = in_scope_probes(scope)
    avail = available_probes(sr=sr)

    per_probe = []
    missing = []
    for pid in probe_ids:
        if not avail.get(pid, False):
            per_probe.append({"probe": pid, "status": "missing"})
            missing.append(pid)
            continue
        try:
            probe = get_probe(pid, sr=sr)
        except MissingProbeSample:
            per_probe.append({"probe": pid, "status": "missing"})
            missing.append(pid)
            continue
        v = run_probe_gate_native(corners_poles, probe, pid, sr)
        v["status"] = "passed" if v["passed"] else "failed"
        per_probe.append(v)

    failed = [p for p in per_probe if p.get("status") == "failed"]
    native_admitted = (not missing) and (not failed)
    return {
        "schema": "trench.gate.cube_verdict_native.v1",
        "interpolation_domain": "pole_freq_radius_linear",
        "numerator_policy": "passthrough_zeros_b0_1",
        "cube_id": f"morph_cube_{cube_decoded['index']:03d}",
        "morpheus_cube_index": cube_decoded["index"],
        "scope": scope,
        "native_edge_pass": native_admitted,
        "missing_samples": missing,
        "per_probe": per_probe,
    }


# -----------------------------------------------------------------------------
# Batch runner
# -----------------------------------------------------------------------------

def load_features(path: Path) -> dict[str, dict]:
    if not path.exists():
        return {}
    with path.open(encoding="utf-8") as fh:
        data = json.load(fh)
    return {e["id"]: e for e in data.get("entries", [])}


def summarize(verdict: dict) -> dict:
    fail_counts: dict[str, int] = {}
    edges_passed = 0
    edges_total = 0
    for p in verdict["per_probe"]:
        if p.get("status") not in ("passed", "failed"):
            continue
        edges_passed += p["edges_passed"]
        edges_total += p["edge_count"]
        for e in p["edges"]:
            for f in e["fails"]:
                fail_counts[f] = fail_counts.get(f, 0) + 1
    return {"edges_passed": edges_passed, "edges_total": edges_total,
            "fail_counts": fail_counts}


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--src", type=Path, default=DEFAULT_VAULT)
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--limit", type=int, default=None)
    parser.add_argument("--start", type=int, default=0)
    parser.add_argument("--scope", default=TRIAGE_SCOPE)
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args(argv)

    args.out.mkdir(parents=True, exist_ok=True)
    with args.src.open(encoding="utf-8") as fh:
        cubes = json.load(fh)["cubes"]
    end = len(cubes) if args.limit is None else min(len(cubes), args.start + args.limit)
    cubes = cubes[args.start:end]
    features = load_features(args.manifest)

    csv_path = args.out / "summary_native.csv"
    fieldnames = ["cube_id", "native_edge_pass", "edges_passed", "edges_total",
                   "fail_dead", "fail_blow", "fail_collapse", "fail_samey",
                   "anomaly_score", "centroid_span_hz", "peak_dynamic_range_db",
                   "rms_dynamic_range_db", "terrain_entropy", "elapsed_s"]
    csv_mode = "w" if (args.force or not csv_path.exists()) else "a"
    with csv_path.open(csv_mode, newline="", encoding="utf-8") as csvfh:
        writer = csv.DictWriter(csvfh, fieldnames=fieldnames)
        if csv_mode == "w":
            writer.writeheader()
        for i, cube in enumerate(cubes):
            cube_id = f"morph_cube_{cube['index']:03d}"
            vpath = args.out / f"{cube_id}.verdict_native.json"
            if vpath.exists() and not args.force:
                print(f"[{i+1}/{len(cubes)}] skip {cube_id}")
                continue
            t0 = time.time()
            v = run_native_edge_gate(cube, scope=args.scope)
            elapsed = time.time() - t0
            vpath.write_text(json.dumps(v, indent=2) + "\n", encoding="utf-8")
            s = summarize(v)
            f = features.get(cube_id, {})
            writer.writerow({
                "cube_id": cube_id,
                "native_edge_pass": v["native_edge_pass"],
                "edges_passed": s["edges_passed"],
                "edges_total": s["edges_total"],
                "fail_dead": s["fail_counts"].get("dead_zone", 0),
                "fail_blow": s["fail_counts"].get("blow_up", 0),
                "fail_collapse": s["fail_counts"].get("abrupt_collapse", 0),
                "fail_samey": s["fail_counts"].get("samey", 0),
                "anomaly_score": f.get("anomaly_score"),
                "centroid_span_hz": f.get("centroid_span_hz"),
                "peak_dynamic_range_db": f.get("peak_dynamic_range_db"),
                "rms_dynamic_range_db": f.get("rms_dynamic_range_db"),
                "terrain_entropy": f.get("terrain_entropy"),
                "elapsed_s": round(elapsed, 2),
            })
            csvfh.flush()
            print(f"[{i+1}/{len(cubes)}] {cube_id} native_pass={v['native_edge_pass']} "
                   f"edges={s['edges_passed']}/{s['edges_total']} fails={s['fail_counts']} "
                   f"{elapsed:.1f}s")
    print(f"done. csv → {csv_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
