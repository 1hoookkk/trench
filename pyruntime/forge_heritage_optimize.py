"""Heritage-space pymoo optimizer v2 — morph-only search, pressure-derived Q.

Search space: 36 params (2 morph endpoints × 6 stages × 3 heritage ints)
Plus 6 discrete pressure selectors for Q derivation.
Total: 42 continuous parameters.

Q axis is derived by applying a gain_shift to the heritage compiler,
matching how the E-MU firmware actually implemented resonance control:
higher Q = larger gain_offset = wider bandwidth separation = sharper peak.

Shipping invariants from SHIPPING.md are hard constraints.
Objectives: talkingness, ridge prominence, morph distance.

Usage:
    python -m pyruntime.forge_heritage_optimize --body unconstrained --generations 100 --pop 80
    python -m pyruntime.forge_heritage_optimize --all --generations 150 --pop 100
"""
from __future__ import annotations

import argparse
import json
import math
import os
from datetime import datetime
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.core.problem import Problem
from pymoo.optimize import minimize
from pymoo.operators.crossover.sbx import SBX
from pymoo.operators.mutation.pm import PM
from pymoo.operators.sampling.rnd import FloatRandomSampling

from pyruntime.body import Body
from pyruntime.constants import SR, NUM_BODY_STAGES
from pyruntime.corner import CornerArray, CornerName, CornerState
from pyruntime.designer_compile import PASSTHROUGH_ENC, _compile_packed
from pyruntime.encode import EncodedCoeffs
from pyruntime.freq_response import cascade_response_db, freq_points
from pyruntime.stage_params import StageParams
from pyruntime.forge_optimize import (
    _talkingness_single, _extract_formants,
    _trajectory_score, _continuity_score,
)

VAULT_DIR = Path(__file__).parent.parent / "vault"
FREQS = freq_points(sr=SR)
TASTE_MODEL_PATH = Path(__file__).parent.parent / "taste_model.json"

_taste_scorer = None
try:
    from pyruntime.forge_optimize import TasteScorer, _extract_taste_features
    if TASTE_MODEL_PATH.is_file():
        _taste_scorer = TasteScorer(TASTE_MODEL_PATH)
except Exception:
    pass


# ---------------------------------------------------------------------------
# Parameter layout
# ---------------------------------------------------------------------------
# Per morph endpoint: 6 stages × 3 params (type, freq, gain) = 18 floats
# 2 endpoints = 36 floats
# Plus 6 Q-intensity params (one per stage, controls gain_shift for Q=1 corners)
# Total: 42 continuous params

N_STAGES = 6
N_PER_STAGE = 3
N_PER_ENDPOINT = N_STAGES * N_PER_STAGE  # 18
N_Q_PARAMS = N_STAGES                     # 6
CROSS_STAGES = [0, 1, 2]                   # formant-driving stages
N_CROSS = len(CROSS_STAGES) * 2             # c1 + c3 delta per stage = 6
N_TOTAL = 2 * N_PER_ENDPOINT + N_Q_PARAMS + N_CROSS  # 48


def _compile_stage(type_id: int, freq: int, gain: int, q_shift: int = 0) -> tuple[StageParams, EncodedCoeffs]:
    """Compile one heritage stage with optional Q shift applied as global_shift."""
    t = max(0, min(3, type_id))
    f = max(0, min(127, freq))
    g = max(0, min(127, gain))
    if t == 0:
        return StageParams.passthrough(), PASSTHROUGH_ENC
    return _compile_packed(t, f, g, shift=q_shift)


def _params_to_body(params: np.ndarray) -> Body:
    """Convert 42-float vector to Body with pressure-derived Q."""
    p_m0 = params[:N_PER_ENDPOINT]          # morph=0 endpoint
    p_m1 = params[N_PER_ENDPOINT:2*N_PER_ENDPOINT]  # morph=1 endpoint
    q_intensities = params[2*N_PER_ENDPOINT:]  # per-stage Q intensity

    corners = []
    for morph, q in [(0.0, 0.0), (0.0, 1.0), (1.0, 0.0), (1.0, 1.0)]:
        p = p_m0 if morph < 0.5 else p_m1
        # Q shift: scale the q_intensity by Q position
        # At Q=0: shift=0 (no modification). At Q=1: shift = q_intensity (full effect)
        # q_intensity ranges -20 to +20 (mapped from optimizer param 0-40)
        stages = []
        pre_encoded = []
        for i in range(N_STAGES):
            base = i * N_PER_STAGE
            t = int(p[base]) % 4
            f = int(round(p[base + 1]))
            g = int(round(p[base + 2]))
            q_shift = int(round((q_intensities[i] - 20.0) * q))
            sp, enc = _compile_stage(t, f, g, q_shift)
            stages.append(sp)
            pre_encoded.append(enc)

        while len(stages) < NUM_BODY_STAGES:
            stages.append(StageParams.passthrough())
            pre_encoded.append(PASSTHROUGH_ENC)

        corners.append(CornerState(stages=stages, boost=4.0, _pre_encoded=pre_encoded))

    ca = CornerArray(a=corners[0], b=corners[1], c=corners[2], d=corners[3])
    return Body(name="opt", corners=ca, boost=4.0)


def _params_to_body_cross(params: np.ndarray) -> Body:
    """48-float vector → Body with Δmq cross-term on stages 0-2.

    The runtime basis is exactly {1, q, m, mq} via bilinear interpolation.
    Adding Δmq to the D corner (M1_Q1) makes the surface non-additive,
    which 99.7% of Morpheus cubes require.

    The shipped body includes the modified D corner, so runtime matches
    what the optimizer evaluated. No custom interpolation needed.

    Params [42:48] = (c3_delta, c4_delta) per CROSS_STAGE.
    c3 = pole frequency (-2*cos(theta)), c4 = pole radius (r²).
    """
    body = _params_to_body(params[:42])
    cross_raw = params[42:48]

    # Modify D corner (M1_Q1) pre-encoded coefficients
    d_corner = body.corners._corners[3]
    d_enc = d_corner.encode()
    new_enc = list(d_enc)
    for idx, si in enumerate(CROSS_STAGES):
        if si < len(new_enc):
            e = new_enc[si]
            new_enc[si] = EncodedCoeffs(
                c0=e.c0, c1=e.c1, c2=e.c2,
                c3=e.c3 + cross_raw[idx * 2],
                c4=e.c4 + cross_raw[idx * 2 + 1],
            )
    d_corner._pre_encoded = new_enc
    return body


# ---------------------------------------------------------------------------
# Evaluation — trajectory-first fitness stack
# ---------------------------------------------------------------------------

_OCC = [i / 4.0 for i in range(5)]  # 5×5 grid positions


def _eval_body(body: Body) -> dict:
    """Trajectory-first fitness stack.

    L1  hard gate   — peak ceiling, dead midpoint
    L2  occupancy   — fraction of 5×5 bins with talk > 0
    L3  trajectory  — F1/F2 path quality over Q = 0.5 morph sweep
    L4  continuity  — smooth formant migration
    L5  flavor      — talk · ridge · coherence (secondary)
    Also returns legacy metrics for shipping constraints.

    Cross-term Δmq is baked into the D corner, so standard bilinear
    interpolation produces the correct surface. No special handling.
    """
    mask_vocal = (FREQS >= 200) & (FREQS <= 5000)
    mask_sub = (FREQS >= 20) & (FREQS <= 60)
    mask_void = (FREQS >= 800) & (FREQS <= 1200)
    mask_hf = (FREQS >= 3000) & (FREQS <= 15000)
    mask_root = (FREQS >= 20) & (FREQS <= 200)

    # -- 5×5 response grid (25 evals) ------------------------------------
    responses: dict[tuple[float, float], np.ndarray] = {}
    for mi in _OCC:
        for qi in _OCC:
            enc = body.corners.interpolate(mi, qi)
            responses[(mi, qi)] = cascade_response_db(enc, FREQS, SR)

    # -- L1: hard gates ---------------------------------------------------
    max_peak = max(float(np.max(db)) for db in responses.values())
    db_mid = responses[(0.5, 0.5)]
    mid_talk = _talkingness_single(db_mid)
    gate_fail = max_peak > 50.0 or mid_talk == 0.0

    # -- L2: occupancy ----------------------------------------------------
    talk_grid = {k: _talkingness_single(db) for k, db in responses.items()}
    dead = sum(1 for t in talk_grid.values() if t == 0.0)
    occupancy = (25 - dead) / 25.0

    # -- L3/L4: trajectory + continuity (Q = 0.5 sweep, 5 pts) -----------
    sweep_keys = [(m, 0.5) for m in _OCC]
    formant_tracks = [_extract_formants(responses[k]) for k in sweep_keys]
    trajectory = _trajectory_score(formant_tracks)
    continuity = _continuity_score(formant_tracks)

    # -- L5: flavour ------------------------------------------------------
    ridge_vals = []
    for db in responses.values():
        pk = float(np.max(db))
        mn = float(np.mean(db))
        ridge_vals.append(max(0.0, min(1.0, (pk - mn) / 40.0)))
    ridge = float(np.mean(ridge_vals))

    tvs = list(talk_grid.values())
    talk_max = max(tvs)
    coherence = (1.0 - float(np.std(tvs) / talk_max)) if talk_max > 0 else 0.0
    alive = [t for t in tvs if t > 0]
    talkingness = float(np.mean(alive)) if alive else 0.0
    flavor = 0.4 * talkingness + 0.3 * ridge + 0.3 * coherence

    # -- morph distance ---------------------------------------------------
    morph_dists = []
    for qi in _OCC:
        d = float(np.sqrt(np.mean((responses[(0.0, qi)] - responses[(1.0, qi)]) ** 2)))
        morph_dists.append(d)
    morph_distance = float(np.mean(morph_dists))

    # -- legacy constraint metrics ----------------------------------------
    sub_peaks = {k: float(np.max(db[mask_sub])) for k, db in responses.items()}
    void_scoops = {k: float(np.max(db[mask_hf])) - float(np.mean(db[mask_void]))
                   for k, db in responses.items()}
    root_peaks = {k: float(np.max(db[mask_root])) for k, db in responses.items()}
    formant_counts = {}
    for k, db in responses.items():
        vocal_db = db[mask_vocal]
        mean_db = float(np.mean(db))
        count = 0
        if len(vocal_db) > 2:
            for i in range(1, len(vocal_db) - 1):
                if vocal_db[i] > vocal_db[i - 1] and vocal_db[i] > vocal_db[i + 1]:
                    if vocal_db[i] > mean_db + 6:
                        count += 1
        formant_counts[k] = count

    return {
        "gate_fail": gate_fail,
        "occupancy": occupancy,
        "trajectory": trajectory,
        "continuity": continuity,
        "flavor": flavor,
        "talkingness": talkingness,
        "ridge": ridge,
        "coherence": coherence,
        "morph_distance": morph_distance,
        "worst_peak": max_peak,
        "sub_peaks": sub_peaks,
        "void_scoops": void_scoops,
        "formant_counts": formant_counts,
        "root_peaks": root_peaks,
    }


# ---------------------------------------------------------------------------
# Shipping constraints
# ---------------------------------------------------------------------------

def _c_speaker_knockerz(m):
    ref = m["sub_peaks"].get((0.0, 0.0), -120.0)
    return max(0.0, max(ref - v for v in m["sub_peaks"].values()) - 6.0)

def _c_bismuth_shrapnel(m):
    return max(0.0, 10.0 - min(m["void_scoops"].values()))

def _c_glottal_snare(m):
    return max(0.0, 2.0 - min(m["formant_counts"].values()))

def _c_phase_molt(m):
    ref = m["root_peaks"].get((0.0, 0.0), -120.0)
    return max(0.0, max(ref - v for v in m["root_peaks"].values()) - 6.0)

def _c_stability(m):
    return max(0.0, m["worst_peak"] - 36.0)

BODY_CONSTRAINTS = {
    "speaker_knockerz": [_c_speaker_knockerz, _c_stability],
    "bismuth_shrapnel": [_c_bismuth_shrapnel, _c_stability],
    "glottal_snare": [_c_glottal_snare, _c_stability],
    "phase_molt": [_c_phase_molt, _c_stability],
    "unconstrained": [_c_stability],
}


# ---------------------------------------------------------------------------
# Pymoo Problem
# ---------------------------------------------------------------------------

class HeritageV2Problem(Problem):
    """48-param search: 36 morph + 6 Q-pressure + 6 midpoint bend.

    Objectives: trajectory, continuity, flavor (all maximised).
    Constraints: hard gate, occupancy ≥ 0.4, plus shipping body constraints.
    """

    def __init__(self, body_name: str = "unconstrained"):
        self.body_name = body_name
        self.shipping = BODY_CONSTRAINTS.get(body_name, BODY_CONSTRAINTS["unconstrained"])

        # Bounds: heritage stages + Q intensity + bend (c1/c3 deltas)
        xl_stage = [0.0, 0.0, 0.0]
        xu_stage = [3.99, 127.0, 127.0]
        xl = (xl_stage * N_STAGES * 2
              + [0.0] * N_Q_PARAMS
              + [-0.3] * N_CROSS)
        xu = (xu_stage * N_STAGES * 2
              + [40.0] * N_Q_PARAMS
              + [0.3] * N_CROSS)

        n_constr = 2 + len(self.shipping)   # gate + occupancy + shipping
        super().__init__(
            n_var=N_TOTAL,
            n_obj=3,
            n_ieq_constr=n_constr,
            xl=np.array(xl), xu=np.array(xu),
        )

    def _evaluate(self, X, out, *args, **kwargs):
        F, G = [], []
        for x in X:
            try:
                body = _params_to_body_cross(x)
                m = _eval_body(body)
                F.append([-m["trajectory"], -m["continuity"], -m["flavor"]])
                g = [
                    1.0 if m["gate_fail"] else -1.0,
                    0.4 - m["occupancy"],
                ]
                g.extend(c(m) for c in self.shipping)
                G.append(g)
            except Exception:
                F.append([0.0, 0.0, 0.0])
                G.append([100.0] * (2 + len(self.shipping)))
        out["F"] = np.array(F)
        out["G"] = np.array(G)


# ---------------------------------------------------------------------------
# Plotting
# ---------------------------------------------------------------------------

def _plot_responses(bodies, names, out_path):
    fig, axes = plt.subplots(1, 3, figsize=(15, 4))
    for body, name in zip(bodies, names):
        for ax, (m, q, title) in zip(axes, [
            (0.0, 0.5, "M0 Q50"), (0.5, 0.5, "M50 Q50"), (1.0, 0.5, "M100 Q50")
        ]):
            enc = body.corners.interpolate(m, q)
            db = cascade_response_db(enc, FREQS, SR)
            ax.semilogx(FREQS, db, label=name, alpha=0.6, linewidth=0.8)
            ax.set_title(title)
            ax.set_xlim(20, SR / 2)
            ax.grid(True, alpha=0.3)
    axes[0].set_ylabel("dB")
    axes[0].legend(fontsize=5, ncol=2)
    plt.tight_layout()
    plt.savefig(out_path, dpi=150)
    plt.close()


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def optimize(body_name: str = "unconstrained",
             generations: int = 100,
             pop_size: int = 80,
             seed: int = 42) -> Path:

    print(f"\n{'='*60}")
    print(f"  Heritage v3: {body_name}")
    print(f"  Search: {N_TOTAL} params (36 morph + 6 Q-pressure + 6 bend)")
    print(f"  Objectives: trajectory, continuity, flavor")
    print(f"  Constraints: gate + occupancy + {[c.__name__ for c in BODY_CONSTRAINTS.get(body_name, [])]}")
    print(f"  Gen: {generations}, Pop: {pop_size}")
    print(f"{'='*60}")

    problem = HeritageV2Problem(body_name)
    algorithm = NSGA2(
        pop_size=pop_size,
        sampling=FloatRandomSampling(),
        crossover=SBX(prob=0.9, eta=10),
        mutation=PM(eta=15),
    )

    result = minimize(problem, algorithm, ("n_gen", generations),
                      seed=seed, verbose=True)

    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = VAULT_DIR / f"hv2_{body_name}_{ts}"
    out_dir.mkdir(parents=True, exist_ok=True)

    if result.X is None or len(result.X) == 0:
        print("  No feasible solutions found.")
        return out_dir

    X = result.X if result.X.ndim == 2 else result.X.reshape(1, -1)
    F = -result.F if result.F.ndim == 2 else (-result.F).reshape(1, -1)

    bodies, names, all_metrics = [], [], []
    for i, x in enumerate(X):
        body = _params_to_body_cross(x)
        name = f"{body_name}_{i:03d}"
        body_named = Body(name=name, corners=body.corners, boost=body.boost)
        m = _eval_body(body_named)

        # Save compiled (trimmed to 6 stages, 4-corner only — bend is optimizer-side)
        compiled = json.loads(body_named.to_compiled_json(provenance="heritage_v3"))
        for kf in compiled["keyframes"]:
            kf["stages"] = kf["stages"][:6]
        (out_dir / f"{name}.compiled.json").write_text(json.dumps(compiled, indent=2))

        # Save keyframe
        (out_dir / f"{name}.keyframe.json").write_text(body_named.to_json())

        bodies.append(body_named)
        names.append(name)
        all_metrics.append(m)

    # Taste scoring
    taste_scores = {}
    if _taste_scorer is not None:
        for body_obj, name in zip(bodies, names):
            try:
                tf = _extract_taste_features(body_obj)
                taste_scores[name] = _taste_scorer.score(tf)
            except Exception:
                taste_scores[name] = 0.0

    # Rank by trajectory-first composite
    def _composite(m):
        return 0.5 * m["trajectory"] + 0.3 * m["continuity"] + 0.2 * m["flavor"]

    ranked = sorted(zip(names, all_metrics, bodies), key=lambda x: -_composite(x[1]))

    print(f"\n  {len(bodies)} Pareto survivors → {out_dir}")
    print(f"\n  Top 15:")
    for name, m, _ in ranked[:15]:
        ts_str = f"taste={taste_scores[name]:.3f}" if name in taste_scores else ""
        print(f"    {name}: traj={m['trajectory']:.2f} cont={m['continuity']:.2f} "
              f"flav={m['flavor']:.2f} occ={m['occupancy']:.2f} "
              f"peak={m['worst_peak']:.1f}dB {ts_str}")

    # Plots
    top_n = min(10, len(ranked))
    _plot_responses([b for _, _, b in ranked[:top_n]],
                    [n for n, _, _ in ranked[:top_n]],
                    out_dir / "responses.png")

    # Ranking file
    with open(out_dir / "ranking.txt", "w") as f:
        f.write(f"# Heritage v3 — {body_name}\n# {len(bodies)} survivors\n\n")
        for i, (name, m, _) in enumerate(ranked):
            f.write(f"{i+1:3d}. {name}: traj={m['trajectory']:.2f} cont={m['continuity']:.2f} "
                    f"flav={m['flavor']:.2f} occ={m['occupancy']:.2f} "
                    f"peak={m['worst_peak']:.1f}dB\n")

    return out_dir


def main():
    parser = argparse.ArgumentParser(description="Heritage v2 optimizer")
    parser.add_argument("--body", default="unconstrained",
                        choices=list(BODY_CONSTRAINTS.keys()))
    parser.add_argument("--generations", type=int, default=100)
    parser.add_argument("--pop", type=int, default=80)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--all", action="store_true")
    args = parser.parse_args()

    if args.all:
        for name in BODY_CONSTRAINTS:
            optimize(name, args.generations, args.pop, args.seed)
    else:
        optimize(args.body, args.generations, args.pop, args.seed)


if __name__ == "__main__":
    main()
