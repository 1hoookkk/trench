"""Overnight iteration loop: generate + rank candidates by P2k role vocabulary distance.

Design intent:
 - Use P2k-derived *role vocabulary* as training data (structural existence proofs),
   not as an audio-metric comparison target.
 - Generate candidates from `pyruntime.forge_generator.BLUEPRINTS`.
 - Score by role distance to `datasets/role_vocab/shipping_role_targets_v1.json`.
 - Fail closed on instability (`r >= 1.0`) via `pyruntime.validator.validate`.
 - Render plots using sonic-table landmarks (`docs/sonic_tables/tables.json`) as evidence rails.

Outputs (additive, no deletes):
  vault/_overnight/<run_id>/
    report.md
    manifest.json
    candidates/<body_key>/*.json
    plots/<body_key>/*_{freq,stage}.png
"""

from __future__ import annotations

import argparse
import datetime as _dt
import json
import random
import subprocess
import sys
from dataclasses import asdict
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"

sys.path.insert(0, str(ROOT))

from pyruntime.body import Body
from pyruntime.constants import SR
from pyruntime.corner import CornerName
import pyruntime.forge_shipping as forge_shipping
from pyruntime.forge_generator import BLUEPRINTS, generate_body
from pyruntime.freq_response import cascade_response_db, freq_points, stage_response
from pyruntime.role_vocab import body_signature, role_distance
from pyruntime.validator import validate


BODY_KEYS: dict[str, str] = {
    "speaker_knockerz": "Speaker Knockerz",
    "aluminum_siding": "Aluminum Siding",
    "small_talk": "Small Talk Ah-Ee",
    "cul_de_sac": "Cul-De-Sac",
}

CORNER_LABELS = [
    (CornerName.A, "M0_Q0"),
    (CornerName.B, "M0_Q100"),
    (CornerName.C, "M100_Q0"),
    (CornerName.D, "M100_Q100"),
]

DEFAULT_VOWEL_INTENT: dict[str, tuple[str, str] | None] = {
    # Name includes the intent; enforce a vowel-shift gesture.
    "small_talk": ("ah", "ee"),
    # Others: discover best-fit vowel endpoints.
    "speaker_knockerz": None,
    "aluminum_siding": None,
    "cul_de_sac": None,
}

DEFAULT_VOWEL_TOL_ST = 6.0
DEFAULT_VOWEL_MIN_DELTA_ST = 1.0


def _slug(text: str) -> str:
    return "".join(ch.lower() if ch.isalnum() else "_" for ch in text).strip("_")


def _is_passthrough(enc) -> bool:
    return (
        abs(enc.c0 - 1.0) < 1e-6
        and abs(enc.c1) < 1e-6
        and abs(enc.c2) < 1e-6
        and abs(enc.c3) < 1e-6
        and abs(enc.c4) < 1e-6
    )


def _load_targets(path: Path) -> dict[str, dict]:
    data = json.loads(path.read_text(encoding="utf-8"))
    targets = data.get("targets", {})
    if not isinstance(targets, dict) or not targets:
        raise ValueError(f"Invalid targets payload: {path}")
    return targets


def _load_landmarks(path: Path) -> list[dict]:
    tables = json.loads(path.read_text(encoding="utf-8"))
    entries = tables.get("landmarks", {}).get("entries", [])
    if not isinstance(entries, list):
        return []
    out = []
    for entry in entries:
        if not isinstance(entry, dict):
            continue
        name = entry.get("name")
        freq_hz = entry.get("freq_hz")
        if isinstance(name, str) and isinstance(freq_hz, (int, float)) and float(freq_hz) > 0:
            out.append({"name": name, "freq_hz": float(freq_hz)})
    return out


def _load_vowels(path: Path) -> dict[str, dict]:
    tables = json.loads(path.read_text(encoding="utf-8"))
    vowels = tables.get("vowels", {})
    if not isinstance(vowels, dict):
        return {}
    out: dict[str, dict] = {}
    for key, entry in vowels.items():
        if not isinstance(key, str) or not isinstance(entry, dict):
            continue
        f1 = entry.get("f1")
        f2 = entry.get("f2")
        label = entry.get("label")
        if not isinstance(f1, (int, float)) or not isinstance(f2, (int, float)):
            continue
        out[key] = {
            "label": label if isinstance(label, str) else key,
            "f1": float(f1),
            "f2": float(f2),
        }
    return out


def _vowel_distance_st(f_hz: float, target_hz: float) -> float:
    if not np.isfinite(f_hz) or not np.isfinite(target_hz) or f_hz <= 0.0 or target_hz <= 0.0:
        return 1.0e9
    return float(12.0 * np.log2(f_hz / target_hz))


def _vowel_distance_2d_st(f1_hz: float, f2_hz: float, vowel: dict) -> float:
    d1 = _vowel_distance_st(f1_hz, vowel["f1"])
    d2 = _vowel_distance_st(f2_hz, vowel["f2"])
    return float(np.sqrt(d1 * d1 + d2 * d2))


def _match_vowel(f1_hz: float, f2_hz: float, vowels: dict[str, dict]) -> dict | None:
    if not vowels:
        return None
    best_key: str | None = None
    best_dist = 1.0e9
    for key, v in vowels.items():
        d1 = _vowel_distance_st(f1_hz, v["f1"])
        d2 = _vowel_distance_st(f2_hz, v["f2"])
        d = float(np.sqrt(d1 * d1 + d2 * d2))
        if d < best_dist:
            best_dist = d
            best_key = key
    if best_key is None:
        return None
    return {
        "key": best_key,
        "label": vowels[best_key]["label"],
        "distance_st": best_dist,
    }


def _formant_tracks_for_morph_sweep(body: Body) -> list[list[float]]:
    tracks: list[list[float]] = []
    for m in forge_shipping._OCC_STEPS:
        enc = body.corners.interpolate(m, 0.5)
        db = forge_shipping.cascade_response_db(enc, forge_shipping.FREQS, forge_shipping.SR)
        tracks.append(forge_shipping._extract_formants(db))
    return tracks


def _vowel_transition_metrics(
    *,
    tracks: list[list[float]],
    vowels: dict[str, dict],
    required: tuple[str, str] | None,
) -> dict:
    start_formants = tracks[0] if tracks else []
    end_formants = tracks[-1] if tracks else []

    start = None
    end = None
    if len(start_formants) >= 2:
        start = _match_vowel(start_formants[0], start_formants[1], vowels)
    if len(end_formants) >= 2:
        end = _match_vowel(end_formants[0], end_formants[1], vowels)

    path_keys: list[str] = []
    for ft in tracks:
        if len(ft) < 2:
            continue
        m = _match_vowel(ft[0], ft[1], vowels)
        if m is None:
            continue
        path_keys.append(str(m["key"]))
    path_unique: list[str] = []
    for key in path_keys:
        if not path_unique or path_unique[-1] != key:
            path_unique.append(key)

    ok = False
    required_eval = None
    if start is not None and end is not None:
        if required is not None:
            req_from = vowels.get(required[0]) if vowels else None
            req_to = vowels.get(required[1]) if vowels else None
            if req_from is not None and req_to is not None and len(start_formants) >= 2 and len(end_formants) >= 2:
                from_dist = _vowel_distance_2d_st(start_formants[0], start_formants[1], req_from)
                to_dist = _vowel_distance_2d_st(end_formants[0], end_formants[1], req_to)
                f1_delta_st = _vowel_distance_st(end_formants[0], start_formants[0])
                f2_delta_st = _vowel_distance_st(end_formants[1], start_formants[1])
                direction = (f1_delta_st <= -DEFAULT_VOWEL_MIN_DELTA_ST) and (f2_delta_st >= DEFAULT_VOWEL_MIN_DELTA_ST)
                required_eval = {
                    "from_dist_st": from_dist,
                    "to_dist_st": to_dist,
                    "distance_sum_st": float(from_dist + to_dist),
                    "f1_delta_st": f1_delta_st,
                    "f2_delta_st": f2_delta_st,
                    "tol_st": DEFAULT_VOWEL_TOL_ST,
                    "min_delta_st": DEFAULT_VOWEL_MIN_DELTA_ST,
                }

                # Small Talk: enforce direction (open/back → close/front), not literal token match.
                if required == ("ah", "ee"):
                    start_ok = start["key"] in {"ah", "aw", "uh", "schwa"}
                    end_ok = end["key"] in {"ee", "ih", "eh", "ae"}
                    ok = start_ok and end_ok and (len(path_unique) >= 2)
                else:
                    within = (from_dist <= DEFAULT_VOWEL_TOL_ST) and (to_dist <= DEFAULT_VOWEL_TOL_ST)
                    ok = within and direction and (len(path_unique) >= 2)
        else:
            ok = len(path_unique) >= 2

    dist_sum = None
    if start is not None and end is not None:
        dist_sum = float(start["distance_st"] + end["distance_st"])

    return {
        "required": {"from": required[0], "to": required[1]} if required is not None else None,
        "required_eval": required_eval,
        "start": start,
        "end": end,
        "distance_sum_st": dist_sum,
        "path_unique": path_unique,
        "ok": ok,
    }


def _run_id(seed: int) -> str:
    ts = _dt.datetime.now().strftime("%Y%m%d_%H%M%S")
    return f"{ts}_seed{seed}"


def _unique_dir(path: Path) -> Path:
    if not path.exists():
        return path
    for i in range(1, 1000):
        candidate = path.with_name(f"{path.name}_{i:03d}")
        if not candidate.exists():
            return candidate
    raise RuntimeError(f"Could not find unique path near {path}")


def _plot_frequency(body: Body, out_path: Path, landmarks: list[dict]) -> None:
    freqs = freq_points(n=1024, sr=SR)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9))
    fig.suptitle(f"{body.name} · Frequency Sweeps", fontsize=14, fontweight="bold", color="#d8d0c2")

    morph_positions = [0.0, 0.17, 0.33, 0.5, 0.67, 0.83, 1.0]
    q_positions = [0.0, 0.25, 0.5, 0.75, 1.0]

    def _decorate(ax):
        ax.set_xlim(20, (SR / 2.0) - 1.0)
        ax.set_ylim(-70, 50)
        ax.set_xlabel("Frequency (Hz)")
        ax.set_ylabel("Magnitude (dB)")
        ax.grid(True, which="both", linewidth=0.6, alpha=0.5)
        for lm in landmarks:
            x = lm["freq_hz"]
            if x < 20.0 or x > (SR / 2.0) - 1.0:
                continue
            ax.axvline(x=x, color="#9fd2c1", linewidth=0.7, alpha=0.12)

    ax = axes[0, 0]
    for i, m in enumerate(morph_positions):
        db = cascade_response_db(body.corners.interpolate(m, 0.0), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"M={m:.0%}", color=plt.cm.cool(i / (len(morph_positions) - 1)))
    ax.set_title("Morph sweep @ Q=0%")
    _decorate(ax)

    ax = axes[0, 1]
    for i, m in enumerate(morph_positions):
        db = cascade_response_db(body.corners.interpolate(m, 1.0), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"M={m:.0%}", color=plt.cm.winter(i / (len(morph_positions) - 1)))
    ax.set_title("Morph sweep @ Q=100%")
    _decorate(ax)

    ax = axes[1, 0]
    for i, q in enumerate(q_positions):
        db = cascade_response_db(body.corners.interpolate(0.0, q), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"Q={q:.0%}", color=plt.cm.autumn(i / (len(q_positions) - 1)))
    ax.set_title("Q sweep @ M=0%")
    _decorate(ax)

    ax = axes[1, 1]
    diag = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]
    for i, d in enumerate(diag):
        db = cascade_response_db(body.corners.interpolate(d, d), freqs, SR)
        ax.semilogx(freqs, db, linewidth=1.2, alpha=0.9, label=f"M=Q={d:.0%}", color=plt.cm.spring(i / (len(diag) - 1)))
    ax.set_title("Diagonal sweep @ M=Q")
    _decorate(ax)

    for ax in axes.flatten():
        ax.legend(fontsize=7, framealpha=0.25, loc="upper right")

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def _plot_stage(body: Body, out_path: Path) -> None:
    freqs = freq_points(n=512, sr=SR)
    center_enc = body.corners.interpolate(0.5, 0.5)

    active_indices = [i for i, s in enumerate(center_enc) if not _is_passthrough(s)]
    if not active_indices:
        active_indices = list(range(len(center_enc)))

    fig, axes = plt.subplots(1, 2, figsize=(15, 6))
    fig.suptitle(f"{body.name} · Stage Analysis", fontsize=14, fontweight="bold", color="#d8d0c2")

    ax = axes[0]
    for i, idx in enumerate(active_indices):
        h = stage_response(center_enc[idx], freqs, SR)
        db = 20.0 * np.log10(np.maximum(np.abs(h), 1e-20))
        ax.semilogx(
            freqs,
            db,
            linewidth=1.0,
            alpha=0.85,
            color=plt.cm.plasma(i / max(1, len(active_indices) - 1)),
            label=f"S{idx + 1}",
        )
    cascade_db = cascade_response_db(center_enc, freqs, SR)
    ax.semilogx(freqs, cascade_db, color="#ffffff", linewidth=2.1, alpha=0.95, label="Cascade")
    ax.set_title("Per-stage responses @ M=50%, Q=50%")
    ax.set_xlim(20, (SR / 2.0) - 1.0)
    ax.set_ylim(-45, 45)
    ax.set_xlabel("Frequency (Hz)")
    ax.set_ylabel("Magnitude (dB)")
    ax.grid(True, which="both", linewidth=0.6, alpha=0.5)
    ax.legend(fontsize=7, ncol=2, framealpha=0.25, loc="upper right")

    peak_freq = np.full((len(center_enc), len(CORNER_LABELS)), np.nan, dtype=float)
    peak_db = np.full((len(center_enc), len(CORNER_LABELS)), np.nan, dtype=float)
    for cidx, (cn, _) in enumerate(CORNER_LABELS):
        enc = body.corners.corner(cn).encode()
        for sidx, stage in enumerate(enc):
            h = stage_response(stage, freqs, SR)
            db = 20.0 * np.log10(np.maximum(np.abs(h), 1e-20))
            peak_i = int(np.argmax(db))
            peak_freq[sidx, cidx] = float(freqs[peak_i])
            peak_db[sidx, cidx] = float(db[peak_i])

    ax = axes[1]
    ax.set_title("Stage peak map (corners)")
    ax.set_xlabel("Corner")
    ax.set_ylabel("Stage")
    ax.set_xlim(-0.5, len(CORNER_LABELS) - 0.5)
    ax.set_ylim(len(center_enc) - 0.5, -0.5)

    for sidx in range(len(center_enc)):
        for cidx, (_, label) in enumerate(CORNER_LABELS):
            f = peak_freq[sidx, cidx]
            db = peak_db[sidx, cidx]
            if not np.isfinite(f) or not np.isfinite(db):
                continue
            alpha = 0.25 + 0.7 * (np.clip((db + 20.0) / 40.0, 0.0, 1.0))
            ax.scatter([cidx], [sidx], s=90, color="#9fd2c1", alpha=float(alpha), edgecolors="none")
            ax.text(cidx, sidx, f"{f:.0f}", color="#0a0a14", fontsize=7, ha="center", va="center", alpha=min(0.95, float(alpha) + 0.2))

    ax.set_xticks(range(len(CORNER_LABELS)))
    ax.set_xticklabels([label for _, label in CORNER_LABELS], rotation=0)
    ax.set_yticks(range(len(center_enc)))
    ax.set_yticklabels([f"S{i + 1}" for i in range(len(center_enc))])
    ax.grid(True, which="both", linewidth=0.5, alpha=0.25)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=170, bbox_inches="tight")
    plt.close(fig)


def _choose_candidates(
    *,
    public_name: str,
    body_key: str,
    target_sig: dict,
    pool_size: int,
    seed: int,
    blueprints: list[str],
    k: int,
    out_dir: Path,
    landmarks: list[dict],
    vowels: dict[str, dict],
) -> list[dict]:
    rng = random.Random(seed)

    scored: list[dict] = []
    attempts = 0
    max_attempts = max(pool_size * 40, pool_size + 128)

    while len(scored) < pool_size and attempts < max_attempts:
        attempts += 1
        archetype = rng.choice(blueprints)
        body_seed = rng.randint(0, 2**31 - 1)
        body = generate_body(archetype, seed=body_seed, solve=True)

        issues = validate(body.corners)
        has_error = any(issue.severity == "error" for issue in issues)
        if has_error:
            continue

        surface = forge_shipping.eval_shipping_body(body, body_key)
        tracks = _formant_tracks_for_morph_sweep(body)
        vowel_required = DEFAULT_VOWEL_INTENT.get(body_key)
        vowel_metrics = _vowel_transition_metrics(tracks=tracks, vowels=vowels, required=vowel_required)

        sig = body_signature(body)
        dist = role_distance(sig, target_sig)

        scored.append(
            {
                "archetype": archetype,
                "seed": body_seed,
                "surface": surface,
                "tracks": tracks,
                "vowels": vowel_metrics,
                "role_distance": dist,
                "issues": [asdict(issue) for issue in issues],
                "body": body,
            }
        )

    def _sort_key(row: dict) -> tuple:
        surface = row["surface"]
        dist = row["role_distance"]
        return (
            0 if not surface.get("gate_fail", True) else 1,
            -float(surface.get("talkingness", 0.0)),
            -float(surface.get("trajectory", 0.0)),
            -float(surface.get("continuity", 0.0)),
            float(dist["composite_distance"]),
        )

    scored.sort(key=_sort_key)
    winners = scored[:k]

    out: list[dict] = []
    for idx, row in enumerate(winners, start=1):
        archetype = row["archetype"]
        body_seed = int(row["seed"])
        dist = row["role_distance"]
        issues = row["issues"]
        surface = row["surface"]
        vowel_metrics = row["vowels"]

        body: Body = row["body"]
        candidate_name = f"{public_name} O/N {out_dir.name} C{idx} ({archetype})"
        body = Body(name=candidate_name, corners=body.corners, boost=body.boost, preset_schema=body.preset_schema)

        candidates_dir = out_dir / "candidates" / body_key
        candidates_dir.mkdir(parents=True, exist_ok=True)
        stem = f"{body_key}__{out_dir.name}__cand_{idx:02d}__{_slug(archetype)}"
        json_path = candidates_dir / f"{stem}.json"
        json_path.write_text(body.to_compiled_json(provenance="overnight-role-loop"), encoding="utf-8")

        plots_dir = out_dir / "plots" / body_key
        plots_dir.mkdir(parents=True, exist_ok=True)
        freq_plot = plots_dir / f"{stem}_freq.png"
        stage_plot = plots_dir / f"{stem}_stage.png"
        _plot_frequency(body, freq_plot, landmarks)
        _plot_stage(body, stage_plot)

        out.append(
            {
                "candidate": f"cand_{idx:02d}",
                "name": body.name,
                "body_key": body_key,
                "public_name": public_name,
                "archetype": archetype,
                "seed": body_seed,
                "path": str(json_path.resolve()),
                "plots": {
                    "frequency": str(freq_plot.resolve()),
                    "stage": str(stage_plot.resolve()),
                },
                "metrics": {
                    "role_distance": dist,
                    "surface": surface,
                    "vowel_transition": vowel_metrics,
                    "validator": {
                        "issues": issues,
                    },
                },
            }
        )

    return out


def main() -> None:
    parser = argparse.ArgumentParser(description="Overnight role-vocab shipping loop (additive outputs).")
    parser.add_argument("--seed", type=int, default=None, help="Optional fixed RNG seed (default: time-based).")
    parser.add_argument("--pool", type=int, default=96, help="Valid candidates generated per shipping body before selecting top-K.")
    parser.add_argument("--topk", type=int, default=4, help="Winners per shipping body.")
    parser.add_argument("--targets", type=Path, default=ROOT / "datasets" / "role_vocab" / "shipping_role_targets_v1.json")
    parser.add_argument("--tables", type=Path, default=ROOT / "docs" / "sonic_tables" / "tables.json")
    parser.add_argument("--out-root", type=Path, default=ROOT / "vault" / "_overnight")
    parser.add_argument("--min-talking", type=float, default=0.0)
    parser.add_argument("--min-trajectory", type=float, default=0.0)
    parser.add_argument("--min-continuity", type=float, default=0.0)
    parser.add_argument(
        "--label-phonemes",
        action="store_true",
        help="After generating the run, label each candidate's corners with nearest PH_### (requires vault/_phonemes).",
    )
    args = parser.parse_args()

    seed = int(args.seed) if args.seed is not None else (int(_dt.datetime.now().timestamp()) & 0x7FFFFFFF)
    if args.pool < 4:
        raise ValueError("--pool must be >= 4")
    if args.topk < 1:
        raise ValueError("--topk must be >= 1")

    targets = _load_targets(args.targets)
    landmarks = _load_landmarks(args.tables)
    vowels = _load_vowels(args.tables)

    blueprints = sorted(BLUEPRINTS.keys())
    if not blueprints:
        raise RuntimeError("No BLUEPRINTS found in pyruntime.forge_generator")

    run_dir = _unique_dir(args.out_root / _run_id(seed))
    run_dir.mkdir(parents=True, exist_ok=True)

    manifest: dict[str, list[dict]] = {}
    summary_rows: list[dict] = []

    for i, (body_key, public_name) in enumerate(BODY_KEYS.items()):
        target = targets.get(public_name)
        if target is None:
            raise KeyError(f"Missing role target for {public_name} in {args.targets}")
        target_sig = target["signature"]

        body_seed = seed + (i + 1) * 10007
        winners = _choose_candidates(
            public_name=public_name,
            body_key=body_key,
            target_sig=target_sig,
            pool_size=args.pool,
            seed=body_seed,
            blueprints=blueprints,
            k=args.topk,
            out_dir=run_dir,
            landmarks=landmarks,
            vowels=vowels,
        )

        manifest[body_key] = winners
        for row in winners:
            dist = row["metrics"]["role_distance"]
            surface = row["metrics"]["surface"]
            vowel_metrics = row["metrics"]["vowel_transition"]
            summary_rows.append(
                {
                    "body": public_name,
                    "candidate": row["name"],
                    "archetype": row["archetype"],
                    "seed": row["seed"],
                    "role_composite_distance": dist["composite_distance"],
                    "role_mismatch_rate": dist["role_mismatch_rate"],
                    "role_mean_l2": dist["mean_l2"],
                    "talkingness": surface.get("talkingness"),
                    "trajectory": surface.get("trajectory"),
                    "continuity": surface.get("continuity"),
                    "peak_db": surface.get("peak_db"),
                    "gate_fail": surface.get("gate_fail"),
                    "vowel_start": (vowel_metrics.get("start") or {}).get("key") if vowel_metrics else None,
                    "vowel_end": (vowel_metrics.get("end") or {}).get("key") if vowel_metrics else None,
                    "vowel_ok": bool(vowel_metrics.get("ok", False)) if vowel_metrics else False,
                    "path": row["path"],
                }
            )

    (run_dir / "manifest.json").write_text(json.dumps(manifest, indent=2), encoding="utf-8")

    report_lines: list[str] = []
    report_lines.append(f"# Overnight role loop · {run_dir.name}")
    report_lines.append("")
    report_lines.append(f"- seed: {seed}")
    report_lines.append(f"- pool per body: {args.pool}")
    report_lines.append(f"- topk per body: {args.topk}")
    report_lines.append(f"- blueprints: {len(blueprints)}")
    report_lines.append(f"- min_talking: {args.min_talking}")
    report_lines.append(f"- min_trajectory: {args.min_trajectory}")
    report_lines.append(f"- min_continuity: {args.min_continuity}")
    report_lines.append(f"- targets: {args.targets}")
    report_lines.append(f"- tables: {args.tables}")
    report_lines.append("")

    by_body: dict[str, list[dict]] = {}
    for row in summary_rows:
        by_body.setdefault(row["body"], []).append(row)

    for body in BODY_KEYS.values():
        report_lines.append(f"## {body}")
        report_lines.append("")
        rows = by_body.get(body, [])
        rows.sort(key=lambda r: (
            0 if (not r.get("gate_fail", True)) else 1,
            0 if r.get("vowel_ok", False) else 1,
            -(float(r.get("talkingness") or 0.0)),
            -(float(r.get("trajectory") or 0.0)),
            -(float(r.get("continuity") or 0.0)),
            float(r["role_composite_distance"]),
        ))
        for idx, row in enumerate(rows, start=1):
            report_lines.append(
                f"{idx}. `{row['candidate']}` · `{row['archetype']}` · seed `{row['seed']}` · "
                f"talk {float(row.get('talkingness') or 0.0):.3f} · traj {float(row.get('trajectory') or 0.0):.3f} · cont {float(row.get('continuity') or 0.0):.3f} · "
                f"vowel {row.get('vowel_start') or '??'}→{row.get('vowel_end') or '??'} ({'OK' if row.get('vowel_ok') else 'NO'}) · "
                f"peak {float(row.get('peak_db') or 0.0):.1f} dB ({'PASS' if not row.get('gate_fail', True) else 'FAIL'}) · "
                f"role {row['role_composite_distance']:.4f} (mismatch {row['role_mismatch_rate']:.3f}, l2 {row['role_mean_l2']:.4f})"
            )
            report_lines.append(f"   - json: `{row['path']}`")
        report_lines.append("")

    (run_dir / "report.md").write_text("\n".join(report_lines).strip() + "\n", encoding="utf-8")

    if args.label_phonemes:
        try:
            subprocess.run(
                [sys.executable, str(TOOLS / "label_overnight_candidates.py"), "--run", str(run_dir)],
                cwd=ROOT,
                check=True,
            )
        except Exception as e:
            (run_dir / "phoneme_labels.ERROR.txt").write_text(str(e), encoding="utf-8")

    print(run_dir)


if __name__ == "__main__":
    main()
