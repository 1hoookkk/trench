"""Edge-only cube gate runner — CUBE_GATE.md §2 + §5 + §6 + §7.

Per the doctrine, the gate renders **probe sources** through the cube (not a
unit impulse). For each in-scope probe, walks the 12 axis-aligned edges,
samples 5 linearly-interpolated positions of the coefficient set, renders
the probe at each, level-normalizes (output RMS == probe RMS), computes
features, and evaluates four fail conditions:

  - dead_zone        : max bin magnitude in any position's output spectrum
                       is more than DEAD_ZONE_DB below the probe's max bin
  - abrupt_collapse  : a single step carries > ABRUPT_COLLAPSE_FRACTION of
                       the total endpoint spectral L2
  - samey            : endpoint spectral L2 < SAMEY_L2_MULT * mean cube-pair
                       L2 for this probe AND endpoint centroid shift below
                       SAMEY_CENTROID_OCTAVE
  - blow_up          : any output bin exceeds probe_max_bin + BLOW_UP_DB

Stage 3 (face coherence), stage 4 (redundancy), and stage 5 source
robustness aggregation across probes are layered on top of this.

Missing probe samples (CUBE_GATE.md §5 admission rule 1) are reported
explicitly and block admission — never silently swapped for another source.
"""
from __future__ import annotations

import argparse
import json
import sys
from itertools import combinations
from pathlib import Path

import numpy as np

from probes import (MissingProbeSample, PROBE_RMS, RUNTIME_SR, available_probes,
                      get_probe, in_scope_probes)
from thumbnail import corner_stage_tuples, f32, run_impulse

CORNER_LABELS = ("c000", "c100", "c010", "c110", "c001", "c101", "c011", "c111")

EDGES = [
    ("c000", "c100"), ("c010", "c110"), ("c001", "c101"), ("c011", "c111"),  # X
    ("c000", "c010"), ("c100", "c110"), ("c001", "c011"), ("c101", "c111"),  # Y
    ("c000", "c001"), ("c100", "c101"), ("c010", "c011"), ("c110", "c111"),  # Z
]
EDGE_POSITIONS = (0.0, 0.25, 0.5, 0.75, 1.0)

# Thresholds — CUBE_GATE.md §2/§9. Refer to the probe's own levels, not to
# absolute dB floors, so real-audio sources are measured in the same frame.
DEAD_ZONE_DB = -24.0                   # max output bin < probe_max_bin + this
BLOW_UP_DB = 24.0                      # any output bin > probe_max_bin + this
ABRUPT_COLLAPSE_FRACTION = 0.60
SAMEY_L2_MULT = 0.10
SAMEY_CENTROID_OCTAVE = 1.0 / 6.0


# -----------------------------------------------------------------------------
# Cascade rendering through a probe (matches trench-core cascade.rs)
# -----------------------------------------------------------------------------

def render_probe(stages: list[tuple[float, float, float, float, float]],
                  boost: float, probe: np.ndarray) -> np.ndarray:
    state = [(0.0, 0.0) for _ in range(len(stages))]
    out = np.empty(probe.size, dtype=np.float32)
    for n in range(probe.size):
        sig = float(probe[n])
        for i, (c0, c1, c2, c3, c4) in enumerate(stages):
            w1, w2 = state[i]
            y = c0 * sig + w1
            new_w1 = c1 * sig - c3 * y + w2
            new_w2 = c2 * sig - c4 * y
            state[i] = (new_w1, new_w2)
            sig = y
        out[n] = f32(sig * boost)
    return out


def level_normalize(output: np.ndarray, probe: np.ndarray) -> np.ndarray:
    """Rescale output RMS to equal probe RMS. Silent output passes through."""
    probe_rms = float(np.sqrt(np.mean(probe.astype(np.float64) ** 2)))
    out_rms = float(np.sqrt(np.mean(output.astype(np.float64) ** 2)))
    if out_rms <= 1e-12:
        return output
    return (output * (probe_rms / out_rms)).astype(np.float32)


def interp_stages(a: list[tuple[float, ...]], b: list[tuple[float, ...]], t: float
                   ) -> list[tuple[float, float, float, float, float]]:
    return [tuple((1.0 - t) * av + t * bv for av, bv in zip(sa, sb))
             for sa, sb in zip(a, b)]


# -----------------------------------------------------------------------------
# Feature extraction
# -----------------------------------------------------------------------------

def magnitude_db(x: np.ndarray, sr: int
                  ) -> tuple[np.ndarray, np.ndarray]:
    spec = np.fft.rfft(x)
    mag = np.maximum(np.abs(spec), 1e-12)
    db = 20.0 * np.log10(mag)
    freqs = np.fft.rfftfreq(x.size, d=1.0 / sr)
    return freqs, db


def spectral_centroid(freqs: np.ndarray, db: np.ndarray) -> float:
    mag = np.power(10.0, db / 20.0)
    weight = mag * mag
    total = weight.sum()
    if total <= 0:
        return 0.0
    return float((freqs * weight).sum() / total)


def l2_spectrum(a_db: np.ndarray, b_db: np.ndarray) -> float:
    return float(np.sqrt(np.mean((a_db - b_db) ** 2)))


def octave_shift(a: float, b: float) -> float:
    if a <= 0 or b <= 0:
        return 0.0
    return abs(np.log2(b / a))


# -----------------------------------------------------------------------------
# Per-probe gate
# -----------------------------------------------------------------------------

def render_cube_corners(cube: dict, probe: np.ndarray, sr: int
                         ) -> dict[str, np.ndarray]:
    """For each corner, render the probe, level-normalize, return spectra (dB)."""
    spectra = {}
    for label in CORNER_LABELS:
        corner = cube["corners"][label]
        stages = corner_stage_tuples(corner)
        boost = float(corner.get("boost", 1.0))
        out = render_probe(stages, boost, probe)
        out = level_normalize(out, probe)
        _, db = magnitude_db(out, sr)
        spectra[label] = db
    return spectra


def mean_pair_l2(spectra: dict[str, np.ndarray]) -> float:
    dists = [l2_spectrum(spectra[a], spectra[b])
              for a, b in combinations(spectra.keys(), 2)]
    return float(np.mean(dists)) if dists else 0.0


def evaluate_edge_under_probe(cube: dict, a_label: str, b_label: str,
                                probe: np.ndarray, probe_max_bin_db: float,
                                mean_pair: float, sr: int) -> dict:
    a_corner = cube["corners"][a_label]
    b_corner = cube["corners"][b_label]
    a_stages = corner_stage_tuples(a_corner)
    b_stages = corner_stage_tuples(b_corner)
    boost_a = float(a_corner.get("boost", 1.0))
    boost_b = float(b_corner.get("boost", 1.0))

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

    # Conditions
    dead_positions = [v for v in max_bins if v < probe_max_bin_db + DEAD_ZONE_DB]
    blown_positions = [v for v in max_bins if v > probe_max_bin_db + BLOW_UP_DB]

    endpoint_l2 = l2_spectrum(spectra_db[0], spectra_db[-1])
    centroid_a = spectral_centroid(freqs, spectra_db[0])
    centroid_b = spectral_centroid(freqs, spectra_db[-1])
    centroid_oct = octave_shift(centroid_a, centroid_b)

    step_l2s = [l2_spectrum(spectra_db[i], spectra_db[i + 1]) for i in range(4)]
    max_step = max(step_l2s) if step_l2s else 0.0
    collapse_frac = (max_step / endpoint_l2) if endpoint_l2 > 1e-9 else 0.0

    samey_l2_hit = endpoint_l2 < SAMEY_L2_MULT * mean_pair
    samey_cent_hit = centroid_oct < SAMEY_CENTROID_OCTAVE
    samey = samey_l2_hit and samey_cent_hit

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


def run_probe_gate(cube: dict, probe: np.ndarray, probe_id: str, sr: int) -> dict:
    probe_max_bin_db = float(magnitude_db(probe, sr)[1].max())
    spectra = render_cube_corners(cube, probe, sr)
    mean_pair = mean_pair_l2(spectra)
    edge_results = [evaluate_edge_under_probe(cube, a, b, probe, probe_max_bin_db,
                                                 mean_pair, sr)
                     for a, b in EDGES]
    all_pass = all(r["passed"] for r in edge_results)
    return {
        "probe": probe_id,
        "passed": all_pass,
        "edge_count": len(edge_results),
        "edges_passed": sum(1 for r in edge_results if r["passed"]),
        "probe_max_bin_db": round(probe_max_bin_db, 2),
        "mean_cube_pair_l2": round(mean_pair, 4),
        "samey_l2_threshold": round(SAMEY_L2_MULT * mean_pair, 4),
        "edges": edge_results,
    }


# -----------------------------------------------------------------------------
# Top-level gate runner
# -----------------------------------------------------------------------------

def run_edge_gate(cube: dict, *, sr: int = RUNTIME_SR,
                    scope_override: str | None = None) -> dict:
    scope = scope_override or cube.get("derived_from", {}).get("scope") \
         or cube.get("scope") or "broadband"
    probe_ids = in_scope_probes(scope)
    avail = available_probes(sr=sr)

    per_probe = []
    missing = []
    for pid in probe_ids:
        if not avail.get(pid, False):
            per_probe.append({"probe": pid, "status": "missing",
                               "reason": "sample not curated — run blocked per "
                                          "CUBE_GATE.md §5 rule 1"})
            missing.append(pid)
            continue
        try:
            probe = get_probe(pid, sr=sr)
        except MissingProbeSample as e:
            per_probe.append({"probe": pid, "status": "missing",
                               "reason": str(e)})
            missing.append(pid)
            continue
        verdict = run_probe_gate(cube, probe, pid, sr)
        verdict["status"] = "passed" if verdict["passed"] else "failed"
        per_probe.append(verdict)

    passed = [p for p in per_probe if p.get("status") == "passed"]
    failed = [p for p in per_probe if p.get("status") == "failed"]
    admitted = (not missing) and (not failed) and len(passed) == len(probe_ids)

    return {
        "schema": "trench.gate.cube_verdict.v1",
        "stage_run": 2,
        "cube_id": cube.get("derived_from", {}).get("id") or cube.get("name", "?"),
        "scope": scope,
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
        "notes": ["stage 3 (faces) + stage 4 (redundancy) not yet implemented; "
                   "stage 5 admission is 'all in-scope probes pass' per §5"],
    }


# -----------------------------------------------------------------------------
# CLI
# -----------------------------------------------------------------------------

def main(argv: list[str]) -> int:
    parser = argparse.ArgumentParser(
        description="Run the edge-only cube gate against a compiled cube surface.",
    )
    parser.add_argument("compiled", help="Path to trench.compiled.cube_surface.v1 JSON.")
    parser.add_argument("--out", default=None, help="Verdict JSON path (stdout if omitted).")
    parser.add_argument("--sr", type=int, default=RUNTIME_SR,
                         help=f"Render sample rate (default: {RUNTIME_SR}).")
    parser.add_argument("--scope", default=None,
                         help="Override the cube's declared scope (for triage).")
    parser.add_argument("--summary", action="store_true",
                         help="Print a human summary to stderr.")
    args = parser.parse_args(argv)

    with open(args.compiled) as f:
        cube = json.load(f)
    if cube.get("schema") != "trench.compiled.cube_surface.v1":
        print(f"error: expected trench.compiled.cube_surface.v1, got "
              f"{cube.get('schema')!r}", file=sys.stderr)
        return 2

    verdict = run_edge_gate(cube, sr=args.sr, scope_override=args.scope)

    text = json.dumps(verdict, indent=2)
    if args.out:
        p = Path(args.out)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text + "\n", encoding="utf-8")
    else:
        sys.stdout.write(text + "\n")

    if args.summary:
        print("\n--- edge gate summary ---", file=sys.stderr)
        print(f"cube: {verdict['cube_id']} (scope={verdict['scope']})", file=sys.stderr)
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
                   f"mean_pair_L2={p['mean_cube_pair_l2']:.3f}",
                  file=sys.stderr)
            for r in p["edges"]:
                if not r["passed"]:
                    fails = ",".join(r["fails"])
                    print(f"      FAIL {r['edge']:<12s} L2={r['endpoint_l2']:.3f} "
                           f"cent={r['centroid_octaves']:.3f}oct "
                           f"collapse={r['abrupt_collapse_fraction']:.3f} "
                           f"[{fails}]", file=sys.stderr)
    return 0 if verdict["admitted"] else 1


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
