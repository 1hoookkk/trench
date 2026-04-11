"""Multi-objective optimization of TRENCH bodies via pymoo.

Per-body objectives and constraints derived from the V1 spec.
Heritage compiler produces the coefficients — pymoo searches the
integer grid (type, freq_a, gain_a, freq_b, gain_b) per stage.

Usage:
    python tools/forge_pymoo.py speaker_knockerz --pop 80 --gen 60
    python tools/forge_pymoo.py ear_bender --pop 80 --gen 60
    python tools/forge_pymoo.py nite_mode --pop 80 --gen 60
    python tools/forge_pymoo.py ah_ay_ee --pop 80 --gen 60
"""
import argparse
import json
import math
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from pymoo.core.problem import ElementwiseProblem
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.repair.rounding import RoundingRepair
from pymoo.operators.sampling.rnd import IntegerRandomSampling
from pymoo.optimize import minimize
from pymoo.termination import get_termination

from pyruntime.designer_compile import legacy_sections_to_four_corner, compile_four_corner_to_body
from pyruntime.analysis import midpoint_audit, morph_trajectory_distance, morph_trajectory_distance_normalized, dense_midpoint_audit
from pyruntime.validator import validate as validate_corners
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.constants import SR
from pyruntime.role_vocab import body_signature, role_distance

OUT_BASE = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "vault")
ROLE_TARGETS_DEFAULT = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "datasets",
    "role_vocab",
    "shipping_role_targets_v1.json",
)

FORGE_KEY_TO_PUBLIC = {
    "speaker_knockerz": "Speaker Knockerz",
    "aluminum_siding": "Aluminum Siding",
    "small_talk_ah_ee": "Small Talk Ah-Ee",
    "cul_de_sac": "Cul-De-Sac",
}


# =========================================================================
# Per-body variable bounds: 6 stages × 5 vars = 30 decision variables
# [type, freq_a, gain_a, freq_b, gain_b] per stage
#
# Bounds are (lower, upper) inclusive integers.
# =========================================================================

BODY_BOUNDS = {
    "speaker_knockerz": [
        # S0-S3: sub cluster -> mids
        {"type": (1, 2), "fa": (2, 14),  "ga": (62, 82), "fb": (38, 72),  "gb": (54, 70)},
        {"type": (1, 2), "fa": (4, 18),  "ga": (60, 80), "fb": (42, 76),  "gb": (56, 72)},
        {"type": (1, 2), "fa": (6, 24),  "ga": (58, 76), "fb": (48, 80),  "gb": (56, 70)},
        {"type": (1, 2), "fa": (10, 30), "ga": (56, 72), "fb": (52, 85),  "gb": (54, 68)},
        # S4-S5: presence
        {"type": (1, 3), "fa": (80, 105), "ga": (48, 66), "fb": (85, 112), "gb": (58, 74)},
        {"type": (1, 3), "fa": (90, 115), "ga": (46, 64), "fb": (92, 118), "gb": (56, 72)},
    ],
    "aluminum_siding": [
        # Permanent 1kHz void. All energy in extreme highs.
        # S0-S1: mid scoop — stays scooped at both endpoints
        {"type": (1, 2), "fa": (68, 82), "ga": (38, 54), "fb": (65, 80), "gb": (36, 52)},
        {"type": (1, 2), "fa": (72, 86), "ga": (40, 56), "fb": (70, 84), "gb": (38, 54)},
        # S2-S3: HF shimmer at rest -> metallic tear at open
        {"type": (1, 3), "fa": (105, 120), "ga": (58, 72), "fb": (112, 127), "gb": (66, 82)},
        {"type": (1, 3), "fa": (108, 122), "ga": (56, 70), "fb": (115, 127), "gb": (64, 80)},
        # S4: ultra-HF push (dog whistle territory)
        {"type": (1, 3), "fa": (115, 127), "ga": (52, 66), "fb": (120, 127), "gb": (68, 84)},
        # S5: low anchor — subtle sub presence, stays put
        {"type": (1, 2), "fa": (8, 25),  "ga": (56, 68), "fb": (8, 25),   "gb": (54, 66)},
    ],
    "small_talk_ah_ee": [
        # Vocal cavity morph: open throat -> strangled digital scream
        # F1: 700Hz (Ah) -> drops/chokes (Ee)
        {"type": (1, 2), "fa": (55, 70), "ga": (64, 80), "fb": (35, 55), "gb": (58, 74)},
        # F2: 1100Hz (Ah) -> 2500Hz (Ee) — dramatic rise
        {"type": (1, 2), "fa": (70, 84), "ga": (62, 78), "fb": (88, 105), "gb": (64, 80)},
        # F3: 2500Hz -> 3000Hz presence bite
        {"type": (1, 2), "fa": (86, 100), "ga": (58, 74), "fb": (92, 108), "gb": (60, 76)},
        # Anti-formant / nasal notch
        {"type": (2, 3), "fa": (78, 98), "ga": (44, 60), "fb": (82, 102), "gb": (44, 60)},
        # HF rasp / digital texture
        {"type": (1, 3), "fa": (100, 118), "ga": (50, 66), "fb": (108, 125), "gb": (56, 72)},
        # Sub anchor — always present
        {"type": (1, 2), "fa": (3, 20), "ga": (58, 72), "fb": (3, 20), "gb": (56, 70)},
    ],
    "cul_de_sac": [
        # Tube-to-comb fracture. Starts thick, fractures at midpoint.
        # S0: root hum anchor — permanent, never disappears
        {"type": (1, 2), "fa": (10, 28), "ga": (64, 80), "fb": (10, 28), "gb": (62, 78)},
        # S1: thick tube resonance at rest -> scattered at open
        {"type": (1, 2), "fa": (30, 50), "ga": (62, 78), "fb": (75, 100), "gb": (58, 74)},
        # S2: mid body -> fractures upward
        {"type": (1, 3), "fa": (45, 65), "ga": (60, 76), "fb": (85, 110), "gb": (56, 72)},
        # S3: upper body -> extreme scatter
        {"type": (1, 3), "fa": (55, 78), "ga": (58, 74), "fb": (95, 120), "gb": (54, 70)},
        # S4: presence -> ultra-HF shards
        {"type": (1, 3), "fa": (70, 92), "ga": (56, 72), "fb": (105, 125), "gb": (52, 68)},
        # S5: ceiling suppressor -> opens to air
        {"type": (1, 3), "fa": (90, 115), "ga": (48, 64), "fb": (110, 127), "gb": (56, 72)},
    ],
}

# Per-body objectives (all minimized by pymoo):
#   f0: -morph_rms (maximize morph motion)
#   f1: midpoint_peak (minimize resonance)
#   f2: body-specific quality metric (minimize)

# Per-body constraints (g <= 0 means feasible):
#   g0: midpoint_peak - 35  (must be <= 35 dB)
#   g1: 3.0 - morph_rms     (must have >= 3 dB morph motion)
#   g2: error_count          (must be 0)


class TrenchBodyProblem(ElementwiseProblem):
    def __init__(self, body_name: str, quality_mode: str = "role_vocab", role_targets_path: str | None = ROLE_TARGETS_DEFAULT):
        self.body_name = body_name
        self.quality_mode = quality_mode
        bounds = BODY_BOUNDS[body_name]
        self.stage_bounds = bounds
        self.freqs = freq_points()
        self.role_target_signature = None
        if role_targets_path and os.path.exists(role_targets_path):
            public_name = FORGE_KEY_TO_PUBLIC.get(body_name)
            if public_name is not None:
                try:
                    payload = json.load(open(role_targets_path, "r", encoding="utf-8"))
                    self.role_target_signature = payload.get("targets", {}).get(public_name, {}).get("signature")
                except Exception:
                    self.role_target_signature = None

        # Build variable bounds: 5 vars per stage × 6 stages = 30
        xl, xu = [], []
        for sb in bounds:
            xl.extend([sb["type"][0], sb["fa"][0], sb["ga"][0], sb["fb"][0], sb["gb"][0]])
            xu.extend([sb["type"][1], sb["fa"][1], sb["ga"][1], sb["fb"][1], sb["gb"][1]])

        super().__init__(
            n_var=30,
            n_obj=3,
            n_ieq_constr=3,
            xl=np.array(xl, dtype=float),
            xu=np.array(xu, dtype=float),
            vtype=int,
        )

    def _decode(self, x):
        sections = []
        for i in range(6):
            base = i * 5
            typ = int(round(x[base]))
            fa = int(round(x[base + 1]))
            ga = int(round(x[base + 2]))
            fb = int(round(x[base + 3]))
            gb = int(round(x[base + 4]))
            sections.append((typ, fa, ga, fb, gb))
        return sections

    def _evaluate(self, x, out, *args, **kwargs):
        sections = self._decode(x)
        name = f"{self.body_name}_opt"

        try:
            template = legacy_sections_to_four_corner(name, sections)
            body = compile_four_corner_to_body(template, boost=4.0)

            # Dense audit (21 morph points × 3 Q points = 63 checks)
            audit = dense_midpoint_audit(body, peak_limit_db=35.0, n_morph=21, n_q=3)

            # Level-normalized morph distance (prevents loudness cheating)
            morph_norm = morph_trajectory_distance_normalized(body)
            morph_raw = morph_trajectory_distance(body)

            issues = validate_corners(body.corners)
            errors = sum(1 for i in issues if i.severity == "error")

            mid_peak = audit["worst_peak_db"]
            norm_rms = morph_norm["normalized_rms_db"]
            raw_rms = morph_raw["mean_rms_db"]

            # Body-specific quality metric
            quality = self._quality_metric(body, sections)

            # Objectives (all minimized) — use NORMALIZED morph distance
            out["F"] = [
                -norm_rms,       # maximize normalized morph motion (no loudness cheating)
                mid_peak,        # minimize dense midpoint peak
                -quality,        # maximize quality (body-specific)
            ]

            # Constraints (g <= 0 is feasible)
            out["G"] = [
                mid_peak - 35.0,       # dense midpoint must be <= 35 dB
                2.0 - norm_rms,        # normalized morph must be >= 2 dB
                float(errors),         # no structural errors
            ]

        except Exception:
            out["F"] = [0.0, 999.0, 0.0]
            out["G"] = [999.0, 999.0, 1.0]

    def _legacy_quality_metric(self, body, sections) -> float:
        """Legacy body-specific quality score (higher = better)."""
        freqs = self.freqs

        if self.body_name == "speaker_knockerz":
            # Sub energy at M0 + mid energy at M100
            enc_m0 = body.corners.interpolate(0.0, 0.0)
            enc_m100 = body.corners.interpolate(1.0, 0.0)
            db_m0 = cascade_response_db(enc_m0, freqs, SR)
            db_m100 = cascade_response_db(enc_m100, freqs, SR)
            sub_mask = (freqs >= 30) & (freqs <= 100)
            mid_mask = (freqs >= 300) & (freqs <= 1500)
            sub_energy = float(np.mean(db_m0[sub_mask])) if np.any(sub_mask) else -60
            mid_energy = float(np.mean(db_m100[mid_mask])) if np.any(mid_mask) else -60
            return sub_energy + mid_energy  # want both strong

        elif self.body_name == "aluminum_siding":
            # Permanent 1kHz void + HF energy increasing with morph
            enc_m0 = body.corners.interpolate(0.0, 0.0)
            enc_m100 = body.corners.interpolate(1.0, 0.0)
            db_m0 = cascade_response_db(enc_m0, freqs, SR)
            db_m100 = cascade_response_db(enc_m100, freqs, SR)
            void_mask = (freqs >= 800) & (freqs <= 1200)
            hf_mask = (freqs >= 5000) & (freqs <= 18000)
            # Want deep void at BOTH endpoints + strong HF at M100
            void_m0 = float(np.mean(db_m0[void_mask])) if np.any(void_mask) else 0
            void_m100 = float(np.mean(db_m100[void_mask])) if np.any(void_mask) else 0
            hf_m100 = float(np.mean(db_m100[hf_mask])) if np.any(hf_mask) else -60
            return -void_m0 - void_m100 + hf_m100  # deep void + bright HF

        elif self.body_name == "small_talk_ah_ee":
            # Formant separation: F1 and F2 should be distinct peaks
            enc_m0 = body.corners.interpolate(0.0, 0.0)
            enc_m100 = body.corners.interpolate(1.0, 0.0)
            db_m0 = cascade_response_db(enc_m0, freqs, SR)
            db_m100 = cascade_response_db(enc_m100, freqs, SR)
            f1_mask = (freqs >= 300) & (freqs <= 900)
            f2_mask = (freqs >= 1500) & (freqs <= 3000)
            valley_mask = (freqs >= 900) & (freqs <= 1500)
            f1_m0 = float(np.max(db_m0[f1_mask])) if np.any(f1_mask) else -60
            f2_m100 = float(np.max(db_m100[f2_mask])) if np.any(f2_mask) else -60
            valley = float(np.mean(db_m0[valley_mask])) if np.any(valley_mask) else -60
            return (f1_m0 - valley) + (f2_m100 - valley)

        elif self.body_name == "cul_de_sac":
            # Tube-to-fracture: concentrated energy at M0, scattered at M100
            enc_m0 = body.corners.interpolate(0.0, 0.0)
            enc_m50 = body.corners.interpolate(0.5, 0.0)
            enc_m100 = body.corners.interpolate(1.0, 0.0)
            db_m0 = cascade_response_db(enc_m0, freqs, SR)
            db_m50 = cascade_response_db(enc_m50, freqs, SR)
            db_m100 = cascade_response_db(enc_m100, freqs, SR)
            # Want: narrow peak at M0, deep null at M50, scattered at M100
            peak_m0 = float(np.max(db_m0) - np.mean(db_m0))  # peakiness
            null_depth = float(-np.max(db_m50))  # deeper null = better
            scatter_m100 = float(np.std(db_m100))  # more variation = more comb-like
            return peak_m0 + null_depth * 0.5 + scatter_m100

        return 0.0

    def _role_vocab_quality(self, body) -> float:
        """Role-vocabulary quality (higher is better = lower distance)."""
        if self.role_target_signature is None:
            return 0.0
        sig = body_signature(body)
        dist = role_distance(sig, self.role_target_signature)
        return -float(dist["composite_distance"])

    def _quality_metric(self, body, sections) -> float:
        """Quality objective selected by mode."""
        legacy = self._legacy_quality_metric(body, sections)
        role_vocab = self._role_vocab_quality(body)

        if self.quality_mode == "legacy":
            return legacy
        if self.quality_mode == "hybrid":
            # Role vocabulary dominates; legacy nudges ties.
            return role_vocab + (0.2 * legacy)
        # default: role-vocabulary objective
        return role_vocab


def run_optimization(
    body_name: str,
    pop_size: int,
    n_gen: int,
    quality_mode: str = "role_vocab",
    role_targets_path: str | None = ROLE_TARGETS_DEFAULT,
):
    problem = TrenchBodyProblem(body_name, quality_mode=quality_mode, role_targets_path=role_targets_path)

    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=IntegerRandomSampling(),
        crossover=SBX(prob=0.9, eta=3.0, vtype=float, repair=RoundingRepair()),
        mutation=PM(eta=3.0, vtype=float, repair=RoundingRepair()),
        eliminate_duplicates=True,
    )

    termination = get_termination("n_gen", n_gen)

    print(f"Running NSGA-II for {body_name}: pop={pop_size}, gen={n_gen}")
    print(f"Search space: 30 integer variables (6 stages × 5 params)")
    print(f"Objectives: max morph motion, min midpoint peak, max quality ({quality_mode})")
    print(f"Constraints: midpoint <= 35dB, morph >= 3dB, no structural errors")
    print()

    res = minimize(problem, algorithm, termination, seed=42, verbose=True)

    # Extract Pareto front
    out_dir = os.path.join(OUT_BASE, f"{body_name}_optimized")
    os.makedirs(out_dir, exist_ok=True)

    if res.X is None:
        print("No feasible solutions found.")
        return

    n_solutions = res.X.shape[0] if res.X.ndim > 1 else 1
    X = res.X if res.X.ndim > 1 else res.X.reshape(1, -1)
    F = res.F if res.F.ndim > 1 else res.F.reshape(1, -1)

    print(f"\n{'='*70}")
    print(f"Pareto front: {n_solutions} solutions")
    print()

    results = []
    for i in range(n_solutions):
        sections = problem._decode(X[i])
        name = f"{body_name}_opt_{i:03d}"

        template = legacy_sections_to_four_corner(name, sections)
        body = compile_four_corner_to_body(template, boost=4.0)

        # Full re-audit with both sparse and dense
        audit_sparse = midpoint_audit(body)
        audit_dense = dense_midpoint_audit(body, peak_limit_db=35.0, n_morph=21, n_q=3)
        morph_raw = morph_trajectory_distance(body)
        morph_norm = morph_trajectory_distance_normalized(body)

        raw_rms = morph_raw["mean_rms_db"]
        norm_rms = morph_norm["normalized_rms_db"]
        dense_peak = audit_dense["worst_peak_db"]
        quality = -F[i][2]

        # Trim to 6 stages for plugin compatibility
        compiled = json.loads(body.to_compiled_json(provenance=f"pymoo-{body_name}"))
        for kf in compiled.get("keyframes", []):
            kf["stages"] = kf["stages"][:6]
        compiled["stages"] = 6

        path = os.path.join(out_dir, f"{name}.json")
        with open(path, "w") as f:
            json.dump(compiled, f, indent=2)

        types_str = "".join(str(s[0]) for s in sections)
        print(f"  {name}  raw={raw_rms:.1f}dB  norm={norm_rms:.1f}dB  peak={dense_peak:+.1f}dB  qual={quality:.1f}  types={types_str}")

        results.append({
            "name": name,
            "raw_morph_db": round(raw_rms, 1),
            "norm_morph_db": round(norm_rms, 1),
            "dense_peak_db": round(dense_peak, 1),
            "quality": round(quality, 1),
            "types": types_str,
            "sections": sections,
            "path": path,
        })

    # Bucket into heroes / edge cases / near-miss legends
    heroes = [r for r in results if r["dense_peak_db"] <= 30.0 and r["norm_morph_db"] >= 5.0 and r["quality"] >= 0]
    edge_cases = [r for r in results if r["dense_peak_db"] <= 35.0 and r["norm_morph_db"] >= 3.0 and r not in heroes]
    near_miss = [r for r in results if r not in heroes and r not in edge_cases]

    heroes.sort(key=lambda r: -r["norm_morph_db"])
    edge_cases.sort(key=lambda r: -r["norm_morph_db"])
    near_miss.sort(key=lambda r: -r["norm_morph_db"])

    summary = {
        "body": body_name,
        "pop_size": pop_size,
        "n_gen": n_gen,
        "pareto_size": n_solutions,
        "heroes": heroes[:10],
        "edge_cases": edge_cases[:10],
        "near_miss_legends": near_miss[:10],
        "worst_dense_peak_db": round(max(r["dense_peak_db"] for r in results), 1) if results else None,
        "all_results": results,
    }
    summary_path = os.path.join(out_dir, "_summary.json")
    with open(summary_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"\n{'='*70}")
    print(f"Generated: {n_solutions}  |  Heroes: {len(heroes)}  |  Edge: {len(edge_cases)}  |  Near-miss: {len(near_miss)}")
    print(f"Worst dense peak: {summary['worst_dense_peak_db']} dB")
    print()

    if heroes:
        print("HEROES (peak<=30, norm_morph>=5, qual>=0):")
        for r in heroes[:5]:
            print(f"  {r['name']}  norm={r['norm_morph_db']:.1f}dB  raw={r['raw_morph_db']:.1f}dB  peak={r['dense_peak_db']:+.1f}dB  qual={r['quality']:.1f}")
    if edge_cases:
        print("EDGE CASES (peak<=35, norm_morph>=3):")
        for r in edge_cases[:5]:
            print(f"  {r['name']}  norm={r['norm_morph_db']:.1f}dB  raw={r['raw_morph_db']:.1f}dB  peak={r['dense_peak_db']:+.1f}dB  qual={r['quality']:.1f}")
    if near_miss:
        print("NEAR-MISS LEGENDS:")
        for r in near_miss[:5]:
            print(f"  {r['name']}  norm={r['norm_morph_db']:.1f}dB  raw={r['raw_morph_db']:.1f}dB  peak={r['dense_peak_db']:+.1f}dB  qual={r['quality']:.1f}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("body", choices=list(BODY_BOUNDS.keys()))
    parser.add_argument("--pop", type=int, default=80)
    parser.add_argument("--gen", type=int, default=60)
    parser.add_argument("--quality-mode", choices=["role_vocab", "legacy", "hybrid"], default="role_vocab")
    parser.add_argument("--role-targets", type=str, default=ROLE_TARGETS_DEFAULT)
    args = parser.parse_args()
    run_optimization(
        args.body,
        args.pop,
        args.gen,
        quality_mode=args.quality_mode,
        role_targets_path=args.role_targets,
    )
