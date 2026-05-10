"""Edge-only pill gate — the 2D (compiled-v1) analogue of edge_gate.py.

A pill has 4 keyframes in a bilinear square (morph × Q). That makes 4
axis-aligned edges. Otherwise the per-probe methodology matches
`edge_gate.py` exactly — render probe through each interpolated position,
level-normalize, extract features, evaluate the same four fail conditions
(dead_zone, blow_up, abrupt_collapse, samey).

Admission rule: every in-scope probe must run and pass on every edge
(CUBE_GATE.md §5). Missing probe samples block admission.
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

from edge_gate import (ABRUPT_COLLAPSE_FRACTION, BLOW_UP_DB, DEAD_ZONE_DB,
                        EDGE_POSITIONS, SAMEY_CENTROID_OCTAVE, SAMEY_L2_MULT,
                        interp_stages, l2_spectrum, level_normalize,
                        magnitude_db, octave_shift, render_probe,
                        spectral_centroid)
from probes import (MissingProbeSample, PROBE_RMS, RUNTIME_SR, available_probes,
                      get_probe, in_scope_probes)

PILL_LABELS = ("M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100")
PILL_EDGES = [
    ("M0_Q0",   "M100_Q0"),     # morph axis, low Q
    ("M0_Q100", "M100_Q100"),   # morph axis, high Q
    ("M0_Q0",   "M0_Q100"),     # Q axis, low morph
    ("M100_Q0", "M100_Q100"),   # Q axis, high morph
]


def pill_keyframes_by_label(pill: dict) -> dict[str, dict]:
    out = {}
    for kf in pill.get("keyframes", []):
        label = kf.get("label")
        if label in PILL_LABELS:
            out[label] = kf
    return out


def keyframe_stage_tuples(kf: dict) -> list[tuple[float, float, float, float, float]]:
    out = []
    for s in kf["stages"]:
        out.append((float(s["c0"]), float(s["c1"]), float(s["c2"]),
                    float(s["c3"]), float(s["c4"])))
    return out


def render_pill_corners(pill: dict, probe: np.ndarray, sr: int
                         ) -> dict[str, np.ndarray]:
    kfs = pill_keyframes_by_label(pill)
    spectra = {}
    for label in PILL_LABELS:
        kf = kfs[label]
        stages = keyframe_stage_tuples(kf)
        boost = float(kf.get("boost", 1.0))
        out = render_probe(stages, boost, probe)
        out = level_normalize(out, probe)
        _, db = magnitude_db(out, sr)
        spectra[label] = db
    return spectra


def mean_pair_l2(spectra: dict[str, np.ndarray]) -> float:
    dists = [l2_spectrum(spectra[a], spectra[b])
              for a, b in combinations(spectra.keys(), 2)]
    return float(np.mean(dists)) if dists else 0.0


def evaluate_pill_edge_under_probe(pill: dict, a_label: str, b_label: str,
                                     probe: np.ndarray, probe_max_bin_db: float,
                                     mean_pair: float, sr: int) -> dict:
    kfs = pill_keyframes_by_label(pill)
    a_stages = keyframe_stage_tuples(kfs[a_label])
    b_stages = keyframe_stage_tuples(kfs[b_label])
    boost_a = float(kfs[a_label].get("boost", 1.0))
    boost_b = float(kfs[b_label].get("boost", 1.0))

    freqs = None
    spectra_db = []
    max_bins = []
    for t in EDGE_POSITIONS:
        stages = interp_stages(a_stages, b_stages, t)
        boost = (1.0 - t) * boost_a + t * boost_b
        out = render_probe(stages, boost, probe)
        out = level_normalize(out, probe)
        fr, db = magnitude_db(out, sr)
        freqs = fr
        spectra_db.append(db)
        max_bins.append(float(db.max()))
    spectra_db = np.array(spectra_db)

    dead_positions = [v for v in max_bins if v < probe_max_bin_db + DEAD_ZONE_DB]
    blown_positions = [v for v in max_bins if v > probe_max_bin_db + BLOW_UP_DB]

    endpoint_l2 = l2_spectrum(spectra_db[0], spectra_db[-1])
    centroid_a = spectral_centroid(freqs, spectra_db[0])
    centroid_b = spectral_centroid(freqs, spectra_db[-1])
    centroid_oct = octave_shift(centroid_a, centroid_b)

    step_l2s = [l2_spectrum(spectra_db[i], spectra_db[i + 1]) for i in range(4)]
    max_step = max(step_l2s) if step_l2s else 0.0
    collapse_frac = (max_step / endpoint_l2) if endpoint_l2 > 1e-9 else 0.0

    samey = (endpoint_l2 < SAMEY_L2_MULT * mean_pair
              and centroid_oct < SAMEY_CENTROID_OCTAVE)

    fails = []
    if dead_positions:
        fails.append("dead_zone")
    if blown_positions:
        fails.append("blow_up")
    if collapse_frac > ABRUPT_COLLAPSE_FRACTION:
        fails.append("abrupt_collapse")
    if samey:
        fails.append("samey")

    return {
        "edge": f"{a_label}-{b_label}",
        "passed": not fails,
        "fails": fails,
        "endpoint_l2": round(endpoint_l2, 4),
        "centroid_octaves": round(centroid_oct, 4),
        "abrupt_collapse_fraction": round(collapse_frac, 4),
        "max_bin_db_per_position": [round(v, 2) for v in max_bins],
        "step_l2s": [round(v, 4) for v in step_l2s],
    }


def run_probe_gate(pill: dict, probe: np.ndarray, probe_id: str, sr: int) -> dict:
    probe_max_bin_db = float(magnitude_db(probe, sr)[1].max())
    spectra = render_pill_corners(pill, probe, sr)
    mean_pair = mean_pair_l2(spectra)
    edge_results = [evaluate_pill_edge_under_probe(pill, a, b, probe,
                                                      probe_max_bin_db, mean_pair, sr)
                     for a, b in PILL_EDGES]
    all_pass = all(r["passed"] for r in edge_results)
    return {
        "probe": probe_id,
        "passed": all_pass,
        "edge_count": len(edge_results),
        "edges_passed": sum(1 for r in edge_results if r["passed"]),
        "probe_max_bin_db": round(probe_max_bin_db, 2),
        "mean_pill_pair_l2": round(mean_pair, 4),
        "samey_l2_threshold": round(SAMEY_L2_MULT * mean_pair, 4),
        "edges": edge_results,
    }


def run_pill_gate(pill: dict, *, sr: int = RUNTIME_SR,
                    scope_override: str | None = None) -> dict:
    scope = scope_override or pill.get("provenance", {}).get("scope") \
         or pill.get("scope") or "broadband"
    probe_ids = in_scope_probes(scope)
    avail = available_probes(sr=sr)

    per_probe = []
    missing = []
    for pid in probe_ids:
        if not avail.get(pid, False):
            per_probe.append({"probe": pid, "status": "missing",
                               "reason": "sample not curated — blocks admission "
                                          "(CUBE_GATE.md §5 rule 1)"})
            missing.append(pid)
            continue
        try:
            probe = get_probe(pid, sr=sr)
        except MissingProbeSample as e:
            per_probe.append({"probe": pid, "status": "missing", "reason": str(e)})
            missing.append(pid)
            continue
        verdict = run_probe_gate(pill, probe, pid, sr)
        verdict["status"] = "passed" if verdict["passed"] else "failed"
        per_probe.append(verdict)

    passed = [p for p in per_probe if p.get("status") == "passed"]
    failed = [p for p in per_probe if p.get("status") == "failed"]
    admitted = (not missing) and (not failed) and len(passed) == len(probe_ids)

    return {
        "schema": "trench.gate.pill_verdict.v1",
        "stage_run": 2,
        "pill_id": pill.get("provenance", {}).get("id") or pill.get("name", "?"),
        "scope": scope,
        "dimensions": 2,
        "in_scope_probes": list(probe_ids),
        "missing_samples": missing,
        "admitted": admitted,
        "per_probe": per_probe,
        "thresholds": {
            "dead_zone_db": DEAD_ZONE_DB,
            "blow_up_db": BLOW_UP_DB,
            "abrupt_collapse_fraction": ABRUPT_COLLAPSE_FRACTION,
            "samey_l2_mult": SAMEY_L2_MULT,
            "samey_centroid_octave": SAMEY_CENTROID_OCTAVE,
            "probe_rms_nominal": PROBE_RMS,
        },
        "notes": ["stage 3 (faces) and stage 4 (redundancy) not applicable to "
                   "pills (2D, no face interior to check); stage 5 admission "
                   "is 'all in-scope probes pass all 4 edges'."],
    }


def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Run the edge-only pill gate against a compiled-v1 pill.",
    )
    parser.add_argument("compiled", help="Path to compiled-v1 pill JSON.")
    parser.add_argument("--out", default=None, help="Verdict JSON path (stdout if omitted).")
    parser.add_argument("--sr", type=int, default=RUNTIME_SR,
                         help=f"Render sample rate (default: {RUNTIME_SR}).")
    parser.add_argument("--scope", default=None,
                         help="Override the pill's declared scope (for triage).")
    parser.add_argument("--summary", action="store_true",
                         help="Print a human summary to stderr.")
    args = parser.parse_args(argv)

    with open(args.compiled) as f:
        pill = json.load(f)
    if pill.get("format") != "compiled-v1":
        print(f"error: expected compiled-v1, got format={pill.get('format')!r}",
              file=sys.stderr)
        return 2

    verdict = run_pill_gate(pill, sr=args.sr, scope_override=args.scope)
    text = json.dumps(verdict, indent=2)
    if args.out:
        p = Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text + "\n", encoding="utf-8")
    else:
        sys.stdout.write(text + "\n")

    if args.summary:
        print("\n--- pill gate summary ---", file=sys.stderr)
        print(f"pill: {verdict['pill_id']} (scope={verdict['scope']})", file=sys.stderr)
        print(f"admitted: {verdict['admitted']}", file=sys.stderr)
        if verdict["missing_samples"]:
            print(f"missing samples (blocks admission): "
                   f"{verdict['missing_samples']}", file=sys.stderr)
        for p in verdict["per_probe"]:
            status = p.get("status", "?")
            if status == "missing":
                print(f"  MISSING  {p['probe']}", file=sys.stderr)
                continue
            print(f"  {status.upper():<8s} {p['probe']:<12s} "
                   f"edges={p['edges_passed']}/{p['edge_count']}  "
                   f"mean_pair_L2={p['mean_pill_pair_l2']:.3f}",
                  file=sys.stderr)
            for r in p["edges"]:
                if not r["passed"]:
                    fails = ",".join(r["fails"])
                    print(f"      FAIL {r['edge']:<20s} L2={r['endpoint_l2']:.3f} "
                           f"cent={r['centroid_octaves']:.3f}oct "
                           f"collapse={r['abrupt_collapse_fraction']:.3f} "
                           f"[{fails}]", file=sys.stderr)
    return 0 if verdict["admitted"] else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
