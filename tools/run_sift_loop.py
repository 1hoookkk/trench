"""Headless sift loop: generate candidates via splice/perturb/cross/target strategies,
score with MCP-style stability + corridor analysis, archive winners.

No server required. Runs standalone against pyruntime modules directly.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import math
import random
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.corner import CornerArray, CornerName
from pyruntime.forge_generator import BLUEPRINTS, generate_body
import pyruntime.forge_shipping as forge_shipping
from pyruntime.validator import validate
from pyruntime.freq_response import cascade_response_db
from pyruntime.constants import SR

STRATEGIES = ["splice", "perturb", "cross", "blueprint"]

BODY_KEYS = {
    "speaker_knockerz": "Speaker Knockerz",
    "aluminum_siding": "Aluminum Siding",
    "small_talk": "Small Talk Ah-Ee",
    "cul_de_sac": "Cul-De-Sac",
}

TABLES_PATH = ROOT / "docs" / "sonic_tables" / "tables.json"


def _bark(hz: float) -> float:
    """Zwicker Bark scale (psychoacoustic frequency axis)."""
    f = max(1.0, float(hz))
    return 13.0 * math.atan(0.00076 * f) + 3.5 * math.atan((f / 7500.0) ** 2)


def _load_sonic_tables(path: Path) -> tuple[dict[str, dict], list[dict]]:
    data = json.loads(path.read_text(encoding="utf-8"))
    vowels = data.get("vowels", {})
    landmarks = data.get("landmarks", {}).get("entries", [])

    out_vowels: dict[str, dict] = {}
    if isinstance(vowels, dict):
        for k, v in vowels.items():
            if not isinstance(k, str) or not isinstance(v, dict):
                continue
            f1 = v.get("f1")
            f2 = v.get("f2")
            if not isinstance(f1, (int, float)) or not isinstance(f2, (int, float)):
                continue
            out_vowels[k] = {
                "label": str(v.get("label", k)),
                "f1": float(f1),
                "f2": float(f2),
            }

    out_landmarks: list[dict] = []
    if isinstance(landmarks, list):
        for entry in landmarks:
            if not isinstance(entry, dict):
                continue
            name = entry.get("name")
            freq_hz = entry.get("freq_hz")
            bw_hz = entry.get("bw_hz", 0.0)
            if isinstance(name, str) and isinstance(freq_hz, (int, float)) and float(freq_hz) > 0.0:
                out_landmarks.append(
                    {
                        "name": name,
                        "freq_hz": float(freq_hz),
                        "bw_hz": float(bw_hz) if isinstance(bw_hz, (int, float)) else 0.0,
                    }
                )

    return out_vowels, out_landmarks


def _match_vowel_bark(f1_hz: float, f2_hz: float, vowels: dict[str, dict]) -> dict | None:
    if not vowels or f1_hz <= 0.0 or f2_hz <= 0.0:
        return None
    f1b = _bark(f1_hz)
    f2b = _bark(f2_hz)
    best_key = None
    best_dist = 1.0e9
    for key, v in vowels.items():
        d1 = f1b - _bark(v["f1"])
        d2 = f2b - _bark(v["f2"])
        d = math.sqrt(d1 * d1 + d2 * d2)
        if d < best_dist:
            best_dist = d
            best_key = key
    if best_key is None:
        return None
    # Treat far matches as unknown (prevents misleading labels).
    if best_dist > 1.5:
        return None
    return {"key": best_key, "label": vowels[best_key]["label"], "distance_bark": float(best_dist)}


def _local_extrema_freqs(
    db: np.ndarray,
    freqs: np.ndarray,
    *,
    kind: str,
    threshold_db: float,
    max_count: int,
) -> list[float]:
    if len(db) < 3:
        return []
    mean_db = float(np.mean(db))
    out: list[tuple[float, float]] = []
    if kind == "peak":
        for i in range(1, len(db) - 1):
            if db[i] > db[i - 1] and db[i] > db[i + 1] and db[i] > mean_db + threshold_db:
                out.append((float(db[i]), float(freqs[i])))
        out.sort(key=lambda p: -p[0])
    elif kind == "notch":
        for i in range(1, len(db) - 1):
            if db[i] < db[i - 1] and db[i] < db[i + 1] and db[i] < mean_db - threshold_db:
                out.append((float(db[i]), float(freqs[i])))
        out.sort(key=lambda p: p[0])
    else:
        raise ValueError(f"unknown kind '{kind}'")

    return [f for _, f in out[:max_count]]


def _landmark_hits(*, freqs_hz: list[float], landmarks: list[dict]) -> list[str]:
    hits: list[tuple[float, str]] = []
    for f in freqs_hz:
        fb = _bark(f)
        for lm in landmarks:
            lf = float(lm["freq_hz"])
            bw = float(lm.get("bw_hz", 0.0))
            tol_hz = max(50.0, bw * 0.5) if bw > 0.0 else 120.0
            lo = max(1.0, lf - tol_hz)
            hi = lf + tol_hz
            tol_bark = max(0.25, (_bark(hi) - _bark(lo)) * 0.5)
            d = abs(fb - _bark(lf))
            if d <= tol_bark:
                hits.append((d, f"landmarks.{lm['name']}"))

    hits.sort(key=lambda x: x[0])
    uniq: list[str] = []
    for _, ref in hits:
        if ref not in uniq:
            uniq.append(ref)
    return uniq


def _shape_labels(body: Body, vowels: dict[str, dict], landmarks: list[dict]) -> dict:
    sweep = [i / 10.0 for i in range(11)]
    formant_tracks: list[list[float]] = []
    vowel_path: list[str] = []
    for m in sweep:
        enc = body.corners.interpolate(m, 0.5)
        db = cascade_response_db(enc, forge_shipping.FREQS, SR)
        ft = forge_shipping._extract_formants(db)
        formant_tracks.append(ft)
        if len(ft) >= 2:
            vm = _match_vowel_bark(ft[0], ft[1], vowels)
            if vm is not None:
                key = str(vm["key"])
                if not vowel_path or vowel_path[-1] != key:
                    vowel_path.append(key)

    start = None
    end = None
    mid = None
    if formant_tracks and len(formant_tracks[0]) >= 2:
        start = _match_vowel_bark(formant_tracks[0][0], formant_tracks[0][1], vowels)
    if formant_tracks and len(formant_tracks[-1]) >= 2:
        end = _match_vowel_bark(formant_tracks[-1][0], formant_tracks[-1][1], vowels)
    mid_idx = len(formant_tracks) // 2
    if formant_tracks and len(formant_tracks[mid_idx]) >= 2:
        mid = _match_vowel_bark(formant_tracks[mid_idx][0], formant_tracks[mid_idx][1], vowels)

    enc_mid = body.corners.interpolate(0.5, 0.5)
    freqs = forge_shipping.FREQS
    db_mid = cascade_response_db(enc_mid, freqs, SR)
    peak_freqs = _local_extrema_freqs(db_mid, freqs, kind="peak", threshold_db=6.0, max_count=6)
    notch_freqs = _local_extrema_freqs(db_mid, freqs, kind="notch", threshold_db=6.0, max_count=6)
    lm_refs = _landmark_hits(freqs_hz=peak_freqs + notch_freqs, landmarks=landmarks)

    refs: list[str] = []
    if mid is not None:
        refs.append(f"vowels.{mid['key']}")
    if start is not None and (mid is None or start["key"] != mid["key"]):
        refs.append(f"vowels.{start['key']}")
    if end is not None and (mid is None or end["key"] != mid["key"]):
        refs.append(f"vowels.{end['key']}")
    refs.extend(lm_refs)

    return {
        "vowel_start": start,
        "vowel_mid": mid,
        "vowel_end": end,
        "vowel_path": vowel_path,
        "table_refs": refs[:8],
    }


# Load P2K cracked bodies for splice/cross source material
def _load_cracked() -> list[Body]:
    cracked_dir = ROOT / "vault" / "_cracked"
    if not cracked_dir.exists():
        return []
    bodies = []
    for f in sorted(cracked_dir.glob("*.json")):
        try:
            b = Body.from_json(str(f))
            bodies.append(b)
        except Exception as e:
            print(f"  skip {f.name}: {e}", file=sys.stderr)
    return bodies


def _splice(a: Body, b: Body, rng: random.Random) -> Body:
    """Take corners from two bodies: A gets morph=0, B gets morph=1."""
    ca = a.corners
    cb = b.corners
    spliced = CornerArray(
        a=ca.corner(CornerName.A),
        b=ca.corner(CornerName.B),
        c=cb.corner(CornerName.C),
        d=cb.corner(CornerName.D),
    )
    return Body(name=f"splice_{a.name[:8]}_{b.name[:8]}", corners=spliced, boost=a.boost)


def _cross(a: Body, b: Body, rng: random.Random) -> Body:
    """Stage-level crossover: stage-wise swap of encoded coefficients."""
    from pyruntime.corner import CornerState
    corners = {}
    for cn in [CornerName.A, CornerName.B, CornerName.C, CornerName.D]:
        ac_state = a.corners.corner(cn)
        bc_state = b.corners.corner(cn)
        ac_enc = ac_state.encode()
        bc_enc = bc_state.encode()
        # Take stages 0-2 from A, 3-5 from B, rest from A
        merged_enc = list(ac_enc[:3]) + list(bc_enc[3:6]) + list(ac_enc[6:])
        # Build a CornerState using pre-encoded bypass (same approach as _splice via pre_encoded)
        corners[cn] = CornerState(
            stages=list(ac_state.stages),  # stages are ignored when _pre_encoded is set
            boost=a.boost,
            _pre_encoded=merged_enc,
        )

    return Body(
        name=f"cross_{a.name[:8]}_{b.name[:8]}",
        corners=CornerArray(
            a=corners[CornerName.A],
            b=corners[CornerName.B],
            c=corners[CornerName.C],
            d=corners[CornerName.D],
        ),
        boost=a.boost,
    )


def _score(body: Body, body_key: str = "vocal") -> dict:
    """Score a body by shape-agnostic surface metrics. Role-independent:
    body_key is kept for signature only — eval_shipping_body ignores it.

    Composite rewards diverse shapes, not only vocal/talky ones:
      - trajectory: formant motion quality
      - continuity: smooth sweep
      - ridge: spectral contrast
      - morph_distance: how far the surface travels
      - ruggedness: structural complexity
    """
    issues = validate(body.corners)
    has_error = any(i.severity == "error" for i in issues)
    if has_error:
        return {"valid": False, "reason": "validation_error"}

    surface = forge_shipping.eval_shipping_body(body, body_key)
    gate_fail = surface.get("gate_fail", True)
    talk = float(surface.get("talkingness", 0.0))
    traj = float(surface.get("trajectory", 0.0))
    cont = float(surface.get("continuity", 0.0))
    peak = float(surface.get("peak_db", 999.0))
    ridge = float(surface.get("ridge", 0.0))
    rugg = float(surface.get("ruggedness", 0.0))
    morph_dist = float(surface.get("morph_distance", 0.0))
    occ = float(surface.get("occupancy", 0.0))
    dyn = float(surface.get("dynamic_range", 0.0))

    # Normalize morph_distance to [0,1]: ~40 dB RMS is a big sweep
    morph_norm = min(1.0, morph_dist / 40.0)
    dyn_norm = min(1.0, dyn / 100.0)

    # Shape-coverage composite: reward structure + motion, not just vocals
    composite = (
        occ * 0.15          # must be alive across the surface
        + traj * 0.20       # formants actually move
        + cont * 0.15       # smoothly, not as jumps
        + ridge * 0.15      # has peaks that stand out
        + morph_norm * 0.20 # surface travels far
        + rugg * 0.10       # has structural complexity
        + talk * 0.05       # tiny vocal bonus, not a dominant factor
    )

    return {
        "valid": True,
        "gate_fail": gate_fail,
        "talkingness": talk,
        "trajectory": traj,
        "continuity": cont,
        "ridge": ridge,
        "ruggedness": rugg,
        "morph_distance": morph_dist,
        "occupancy": occ,
        "dynamic_range": dyn,
        "peak_db": peak,
        "composite": composite,
    }


def main():
    parser = argparse.ArgumentParser(description="Headless sift loop.")
    parser.add_argument("--count", type=int, default=256, help="Total candidates to generate.")
    parser.add_argument("--topk", type=int, default=8, help="Winners to keep.")
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--out-root", type=Path, default=ROOT / "vault" / "_diagnostics")
    parser.add_argument("--tables", type=Path, default=TABLES_PATH)
    parser.add_argument("--fill", action="store_true", help="Fill to topk even if it repeats signatures.")
    parser.add_argument("--max-per-vowel", type=int, default=2, help="Max winners per mid vowel key.")
    args = parser.parse_args()

    seed = args.seed if args.seed else (int(_dt.datetime.now().timestamp()) & 0x7FFFFFFF)
    rng = random.Random(seed)

    cracked = _load_cracked()
    bp_keys = sorted(BLUEPRINTS.keys())
    vowels, landmarks = _load_sonic_tables(args.tables)

    if not cracked:
        print("WARNING: no cracked bodies found, splice/cross disabled")

    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    out_dir = args.out_root / f"sift_{ts}_seed{seed}"
    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / "candidates").mkdir(exist_ok=True)

    all_results: list[dict] = []

    for i in range(args.count):
        strategy = rng.choice(STRATEGIES)

        body = None
        try:
            if strategy == "blueprint":
                arch = rng.choice(bp_keys)
                body_seed = rng.randint(0, 2**31 - 1)
                body = generate_body(arch, seed=body_seed, solve=True)
                body = Body(name=f"bp_{arch}_{i:04d}", corners=body.corners, boost=body.boost)

            elif strategy == "splice" and len(cracked) >= 2:
                a, b = rng.sample(cracked, 2)
                body = _splice(a, b, rng)

            elif strategy == "cross" and len(cracked) >= 2:
                a, b = rng.sample(cracked, 2)
                body = _cross(a, b, rng)

            elif strategy == "perturb" and cracked:
                # Not implemented yet — skip to blueprint fallback
                arch = rng.choice(bp_keys)
                body_seed = rng.randint(0, 2**31 - 1)
                body = generate_body(arch, seed=body_seed, solve=True)
                body = Body(name=f"bp_{arch}_{i:04d}", corners=body.corners, boost=body.boost)
                strategy = "blueprint"

            else:
                arch = rng.choice(bp_keys)
                body_seed = rng.randint(0, 2**31 - 1)
                body = generate_body(arch, seed=body_seed, solve=True)
                body = Body(name=f"bp_{arch}_{i:04d}", corners=body.corners, boost=body.boost)
                strategy = "blueprint"

        except Exception as e:
            continue

        if body is None:
            continue

        # Score against first shipping key as baseline
        score = _score(body, "speaker_knockerz")
        if not score.get("valid"):
            continue
        # NOTE: gate_fail only rejects blueprint-forged candidates (they should
        # hit the +36dB target). Splice/cross from cracked P2K can legitimately
        # exceed that — let the composite ranker handle them.
        if strategy == "blueprint" and score.get("gate_fail", True):
            continue

        all_results.append({
            "index": i,
            "name": body.name,
            "strategy": strategy,
            "score": score,
            "body": body,
        })

    # --- Coverage sampler: keep one specimen per distinct shape signature ---
    # Signature = (vowel_mid, vowel_path, centroid_bucket, dynamic_range_bucket,
    #              ridge_bucket, primary_landmark)
    # Buckets are coarse so near-duplicates collapse to the same key.

    def _bucket(v: float, step: float) -> int:
        return int(v / step)

    bank: dict[tuple, dict] = {}

    for r in all_results:
        labels = _shape_labels(r["body"], vowels, landmarks)
        r["labels"] = labels
        s = r["score"]

        mid_key = (labels.get("vowel_mid") or {}).get("key") or "none"
        v_path = "->".join(labels.get("vowel_path", [])) if labels.get("vowel_path") else ""
        primary_ref = (labels.get("table_refs") or ["none"])[0]

        # Shape fingerprint buckets
        morph_b = _bucket(s.get("morph_distance", 0.0), 5.0)     # 5 dB RMS steps
        ridge_b = _bucket(s.get("ridge", 0.0), 0.1)              # 0.1 steps
        rugg_b = _bucket(s.get("ruggedness", 0.0), 0.1)
        dyn_b = _bucket(s.get("dynamic_range", 0.0), 20.0)       # 20 dB steps

        signature = (
            str(mid_key),
            str(v_path),
            str(primary_ref),
            morph_b,
            ridge_b,
            rugg_b,
            dyn_b,
        )

        # Keep the first specimen for each bucket (no ranking — collect, don't rank)
        if signature not in bank:
            bank[signature] = r

    # Output: every unique shape, capped only by --topk if specified
    all_specimens = list(bank.values())

    # Optional cap. For a true shape bank, pass --topk 99999.
    if args.topk > 0 and len(all_specimens) > args.topk:
        # Sort deterministically by signature so the cap is stable
        all_specimens.sort(key=lambda r: r["name"])
        winners = all_specimens[: args.topk]
    else:
        winners = all_specimens

    # Write winners
    lines = [f"# Sift loop · {out_dir.name}", "", f"- seed: {seed}", f"- count: {args.count}",
             f"- topk: {args.topk}", f"- strategies: {STRATEGIES}", f"- cracked sources: {len(cracked)}", ""]

    manifest = {
        "seed": seed,
        "count": args.count,
        "topk": args.topk,
        "tables": str(args.tables),
        "winners": [],
    }

    for idx, w in enumerate(winners, 1):
        s = w["score"]
        labels = w.get("labels") or {}
        v_mid = (labels.get("vowel_mid") or {}).get("key", "none")
        v_path = "->".join(labels.get("vowel_path", [])) if labels.get("vowel_path") else ""
        refs = ", ".join(labels.get("table_refs", [])) if labels.get("table_refs") else ""
        line = (f"{idx}. `{w['name']}` · {w['strategy']} · "
                f"talk {s['talkingness']:.3f} · traj {s['trajectory']:.3f} · cont {s['continuity']:.3f} · "
                f"peak {s['peak_db']:.1f} dB · composite {s['composite']:.3f} · "
                f"vowel {v_mid} · path {v_path}")
        lines.append(line)
        if refs:
            lines.append(f"   - refs: {refs}")

        # Save candidate JSON
        fp = out_dir / "candidates" / f"{w['name']}.json"
        fp.write_text(w["body"].to_compiled_json(provenance="sift-loop"), encoding="utf-8")
        lines.append(f"   - json: `{fp}`")

        manifest["winners"].append({
            "rank": idx,
            "name": w["name"],
            "strategy": w["strategy"],
            "path": str(fp),
            "score": s,
            "labels": labels,
        })

    lines.append("")

    # Strategy distribution
    from collections import Counter
    strat_counts = Counter(r["strategy"] for r in all_results)
    lines.append("## Strategy distribution")
    for k, v in strat_counts.most_common():
        lines.append(f"- {k}: {v}")

    (out_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(out_dir)


if __name__ == "__main__":
    main()
