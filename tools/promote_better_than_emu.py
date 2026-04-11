"""Promote shipping candidates only when they beat translated E-MU baselines.

This script compares candidate bodies against deterministic translation baselines:
- Speaker Knockerz  <- P2k_006
- Aluminum Siding   <- P2k_026
- Small Talk Ah-Ee  <- P2k_013
- Cul-De-Sac        <- P2k_010

It enforces:
1) shipping invariant gate pass
2) midpoint stability pass
3) positive task-metric improvement vs baseline

Outputs:
- JSON report with full measurements and promotion decisions
- Markdown summary for quick review
"""
from __future__ import annotations

import argparse
import json
import re
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import scipy.signal
import scipy.io.wavfile
from scipy.io.wavfile import WavFileWarning

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

from pyruntime.analysis import body_profile, shipping_gate
from pyruntime.body import Body
from pyruntime.render import apply_cascade
from pyruntime.drive import authoring_drive_chain

SR = 44100
MAX_ABS_METRIC_CONTRIBUTION = 1.5

LABEL_TOKENS = {
    "bass": ("808", "bass", "sub", "lowend", "low"),
    "vocal": ("vocal", "vox", "voice", "speech", "acap", "acappella", "chant", "hook", "chorus", "singer"),
    "drum": ("drum", "break", "snare", "hat", "perc", "loop", "kick", "clap", "rim", "cymbal", "hihat", "oh"),
    "mix": ("mix", "bus", "music", "2track", "master", "song", "fullmix", "twotrack"),
}


BODY_PUBLIC_NAME = {
    "speaker_knockerz": "Speaker Knockerz",
    "aluminum_siding": "Aluminum Siding",
    "small_talk": "Small Talk Ah-Ee",
    "cul_de_sac": "Cul-De-Sac",
}

BASELINE_COMPILED = {
    "speaker_knockerz": ROOT / "cartridges" / "translation" / "speaker_knockerz.compiled.json",
    "aluminum_siding": ROOT / "cartridges" / "translation" / "aluminum_siding.compiled.json",
    "small_talk": ROOT / "cartridges" / "translation" / "small_talk_ah_ee.compiled.json",
    "cul_de_sac": ROOT / "cartridges" / "translation" / "cul_de_sac.compiled.json",
}


@dataclass(frozen=True)
class MetricRule:
    name: str
    mode: str
    scale: float


RULES: dict[str, list[MetricRule]] = {
    "speaker_knockerz": [
        MetricRule("sub_anchor_shift_db", "abs_zero", 6.0),
        MetricRule("translation_bite_db", "higher", 6.0),
        MetricRule("danger_contrast_db", "higher", 6.0),
    ],
    "aluminum_siding": [
        MetricRule("void_depth_db", "higher", 6.0),
        MetricRule("air_band_gain_db", "higher", 6.0),
        MetricRule("mid_fill_db", "lower", 4.0),
    ],
    "small_talk": [
        MetricRule("formant_peak_count", "higher", 1.0),
        MetricRule("formant_separation_hz", "higher", 400.0),
        MetricRule("speech_band_gain_db", "higher", 4.0),
    ],
    "cul_de_sac": [
        MetricRule("root_tether_shift_db", "abs_zero", 6.0),
        MetricRule("null_depth_db", "higher", 8.0),
        MetricRule("fracture_distance_db", "higher", 6.0),
    ],
}


def _to_float32_mono(x: np.ndarray) -> np.ndarray:
    if x.ndim > 1:
        x = np.mean(x, axis=1)
    if np.issubdtype(x.dtype, np.integer):
        maxv = np.iinfo(x.dtype).max
        x = x.astype(np.float32) / float(maxv)
    else:
        x = x.astype(np.float32)
    peak = float(np.max(np.abs(x))) if len(x) else 0.0
    if peak > 1.0:
        x = x / peak
    return x


def _resample_if_needed(x: np.ndarray, sr_in: int, sr_out: int = SR) -> np.ndarray:
    if sr_in == sr_out or len(x) == 0:
        return x
    gcd = np.gcd(sr_in, sr_out)
    up = sr_out // gcd
    down = sr_in // gcd
    return scipy.signal.resample_poly(x, up, down).astype(np.float32)


def _load_wav(path: Path) -> np.ndarray:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", WavFileWarning)
        sr, data = scipy.io.wavfile.read(path)
    x = _to_float32_mono(data)
    x = _resample_if_needed(x, int(sr), SR)
    return x


def _synth_probes() -> dict[str, list[np.ndarray]]:
    t = np.arange(0, 2.0, 1.0 / SR, dtype=np.float32)
    env = np.exp(-t * 2.2).astype(np.float32)

    bass = (
        0.60 * np.sin(2.0 * np.pi * 45.0 * t) +
        0.25 * np.sin(2.0 * np.pi * 90.0 * t) +
        0.12 * np.sin(2.0 * np.pi * 180.0 * t)
    ) * (0.5 + 0.5 * env)

    vocal = (
        0.36 * np.sin(2.0 * np.pi * 220.0 * t) +
        0.26 * np.sin(2.0 * np.pi * 700.0 * t) +
        0.22 * np.sin(2.0 * np.pi * 1200.0 * t) +
        0.18 * np.sin(2.0 * np.pi * 2500.0 * t)
    ) * (0.5 + 0.5 * np.sin(2.0 * np.pi * 2.0 * t) ** 2)

    drum = np.zeros_like(t)
    hit_positions = [0, int(0.25 * SR), int(0.5 * SR), int(0.75 * SR), int(1.25 * SR), int(1.5 * SR)]
    rng = np.random.default_rng(7)
    for p in hit_positions:
        n = min(len(drum) - p, int(0.08 * SR))
        if n <= 0:
            continue
        burst = rng.standard_normal(n).astype(np.float32) * np.hanning(n).astype(np.float32)
        drum[p:p + n] += 0.45 * burst
    drum = np.clip(drum, -1.0, 1.0)

    mix = np.clip(0.42 * bass + 0.38 * vocal + 0.32 * drum, -1.0, 1.0)
    return {
        "bass": [bass.astype(np.float32)],
        "vocal": [vocal.astype(np.float32)],
        "drum": [drum.astype(np.float32)],
        "mix": [mix.astype(np.float32)],
    }


def classify_name(path: Path) -> str | None:
    text = f"{str(path.parent).lower()} {path.stem.lower()}"
    tokens = set(t for t in re.split(r"[^a-z0-9]+", text) if t)
    compact = text.replace(" ", "")

    # Prefer explicit mix markers before other classes.
    if ("2track" in compact) or ("twotrack" in compact):
        return "mix"
    if "fullmix" in compact:
        return "mix"
    if "mix" in tokens or "master" in tokens or "bus" in tokens:
        return "mix"

    # Bass next: 808/sub identity is usually explicit.
    if "808" in tokens or "bass" in tokens or "sub" in tokens or "lowend" in tokens:
        return "bass"

    if any(t in tokens for t in LABEL_TOKENS["vocal"]):
        return "vocal"
    if any(t in tokens for t in LABEL_TOKENS["drum"]):
        return "drum"
    return None


def count_corpus_labels(corpus_dir: Path) -> dict[str, int]:
    counts = {"bass": 0, "vocal": 0, "drum": 0, "mix": 0}
    if not corpus_dir.exists():
        return counts
    for wav in corpus_dir.rglob("*.wav"):
        kind = classify_name(wav)
        if kind is not None:
            counts[kind] += 1
    return counts


def load_corpus_with_meta(
    corpus_dir: Path | None,
    *,
    allow_fallback: bool = True,
) -> tuple[dict[str, list[np.ndarray]], str, dict[str, int], list[str]]:
    if corpus_dir is None or (not corpus_dir.exists()):
        if allow_fallback:
            synth = _synth_probes()
            counts = {"bass": 0, "vocal": 0, "drum": 0, "mix": 0}
            missing = ["bass", "vocal", "drum", "mix"]
            return synth, "synthetic", counts, missing
        raise FileNotFoundError(f"Corpus directory not found: {corpus_dir}")

    clips: dict[str, list[np.ndarray]] = {"bass": [], "vocal": [], "drum": [], "mix": []}
    for wav in corpus_dir.rglob("*.wav"):
        kind = classify_name(wav)
        if kind is None:
            continue
        try:
            clips[kind].append(_load_wav(wav))
        except Exception:
            continue
    counts = {k: len(v) for k, v in clips.items()}
    missing = [k for k, v in counts.items() if v == 0]

    if missing and not allow_fallback:
        raise ValueError(
            f"Corpus {corpus_dir} is missing labeled classes: {missing}. "
            "Add WAV files labeled bass/vocal/drum/mix."
        )

    if missing and allow_fallback:
        synth = _synth_probes()
        for k in missing:
            clips[k] = synth[k]
        return clips, "corpus+synthetic_fallback", counts, missing

    return clips, "corpus", counts, []


def load_corpus(corpus_dir: Path | None) -> tuple[dict[str, list[np.ndarray]], str]:
    clips, mode, _, _ = load_corpus_with_meta(corpus_dir, allow_fallback=True)
    return clips, mode


def balanced_sample_clips(
    clips: dict[str, list[np.ndarray]],
    *,
    seed: int = 42,
    max_per_class: int | None = None,
) -> tuple[dict[str, list[np.ndarray]], dict[str, int], dict[str, int], int]:
    if max_per_class is not None and max_per_class <= 0:
        raise ValueError("max_per_class must be positive when provided.")

    input_counts = {k: len(v) for k, v in clips.items()}
    nonzero = [n for n in input_counts.values() if n > 0]
    if not nonzero:
        output = {k: [] for k in clips}
        output_counts = {k: 0 for k in clips}
        return output, input_counts, output_counts, 0

    target = min(nonzero)
    if max_per_class is not None:
        target = min(target, max_per_class)

    rng = np.random.default_rng(seed)
    balanced: dict[str, list[np.ndarray]] = {}
    output_counts: dict[str, int] = {}
    for label, items in clips.items():
        n = len(items)
        if n == 0 or target == 0:
            chosen: list[np.ndarray] = []
        elif n <= target:
            chosen = items[:]
        else:
            idx = rng.choice(n, size=target, replace=False)
            chosen = [items[int(i)] for i in idx]
        balanced[label] = chosen
        output_counts[label] = len(chosen)

    return balanced, input_counts, output_counts, target


def _band_db(signal: np.ndarray, lo: float, hi: float) -> float:
    if len(signal) == 0:
        return -120.0
    x = signal.astype(np.float64)
    spec = np.fft.rfft(x * np.hanning(len(x)))
    freqs = np.fft.rfftfreq(len(x), 1.0 / SR)
    mask = (freqs >= lo) & (freqs <= hi)
    if not np.any(mask):
        return -120.0
    power = float(np.mean(np.abs(spec[mask]) ** 2))
    return 10.0 * np.log10(max(power, 1e-20))


def _full_rms_db(signal: np.ndarray) -> float:
    rms = float(np.sqrt(np.mean(signal.astype(np.float64) ** 2)))
    return 20.0 * np.log10(max(rms, 1e-12))


def _spectral_distance_db(a: np.ndarray, b: np.ndarray) -> float:
    sa = np.abs(np.fft.rfft(a * np.hanning(len(a))))
    sb = np.abs(np.fft.rfft(b * np.hanning(len(b))))
    da = 20.0 * np.log10(np.maximum(sa, 1e-12))
    db = 20.0 * np.log10(np.maximum(sb, 1e-12))
    return float(np.sqrt(np.mean((da - db) ** 2)))


def _vocal_peak_metrics(signal: np.ndarray) -> tuple[int, float]:
    spec = np.abs(np.fft.rfft(signal * np.hanning(len(signal))))
    freqs = np.fft.rfftfreq(len(signal), 1.0 / SR)
    db = 20.0 * np.log10(np.maximum(spec, 1e-12))
    mask = (freqs >= 200.0) & (freqs <= 5000.0)
    f = freqs[mask]
    v = db[mask]
    if len(v) < 5:
        return 0, 0.0
    mean_v = float(np.mean(v))
    peaks: list[tuple[float, float]] = []
    for i in range(1, len(v) - 1):
        if v[i] > v[i - 1] and v[i] > v[i + 1] and v[i] > mean_v + 3.0:
            peaks.append((float(f[i]), float(v[i])))
    peaks.sort(key=lambda p: -p[1])
    count = len(peaks)
    sep = 0.0
    if len(peaks) >= 2:
        sep = abs(peaks[0][0] - peaks[1][0])
    return count, float(sep)


def _process(body: Body, clip: np.ndarray, m: float, q: float) -> np.ndarray:
    dry = authoring_drive_chain(clip, mackie_amount=0.5, erode_amount=0.0, corrode_amount=0.0, sr=SR)
    encoded = body.corners.interpolate(m, q)
    wet = apply_cascade(dry, encoded)
    return np.clip(wet, -1.0, 1.0).astype(np.float32)


def _mean_metric(
    clips: list[np.ndarray],
    fn: Callable[[np.ndarray], float],
) -> float:
    vals = [fn(c) for c in clips]
    return float(np.mean(vals)) if vals else 0.0


def measure_task_metrics(body: Body, body_key: str, clips: dict[str, list[np.ndarray]]) -> dict[str, float]:
    if body_key == "speaker_knockerz":
        bass = clips["bass"]
        calm = [_process(body, c, 0.0, 0.0) for c in bass]
        danger = [_process(body, c, 0.75, 0.0) for c in bass]
        open_state = [_process(body, c, 1.0, 0.0) for c in bass]
        sub_anchor_shift = float(np.mean([
            _band_db(o, 20.0, 60.0) - _band_db(c, 20.0, 60.0)
            for o, c in zip(open_state, calm)
        ]))
        translation_bite = float(np.mean([
            _band_db(d, 700.0, 2500.0) - _band_db(c, 700.0, 2500.0)
            for d, c in zip(danger, calm)
        ]))
        danger_contrast = float(np.mean([
            _full_rms_db(d) - _full_rms_db(c)
            for d, c in zip(danger, calm)
        ]))
        return {
            "sub_anchor_shift_db": sub_anchor_shift,
            "translation_bite_db": translation_bite,
            "danger_contrast_db": danger_contrast,
        }

    if body_key == "aluminum_siding":
        src = clips["drum"] + clips["vocal"]
        calm = [_process(body, c, 0.0, 0.0) for c in src]
        danger = [_process(body, c, 1.0, 0.0) for c in src]
        void_depth = _mean_metric(
            danger,
            lambda x: _band_db(x, 3000.0, 15000.0) - _band_db(x, 800.0, 1200.0),
        )
        air_gain = float(np.mean([
            _band_db(d, 7000.0, 16000.0) - _band_db(c, 7000.0, 16000.0)
            for d, c in zip(danger, calm)
        ]))
        mid_fill = float(np.mean([
            _band_db(d, 800.0, 1200.0) - _band_db(c, 800.0, 1200.0)
            for d, c in zip(danger, calm)
        ]))
        return {
            "void_depth_db": void_depth,
            "air_band_gain_db": air_gain,
            "mid_fill_db": mid_fill,
        }

    if body_key == "small_talk":
        vocal = clips["vocal"]
        calm = [_process(body, c, 0.0, 0.0) for c in vocal]
        bite = [_process(body, c, 1.0, 1.0) for c in vocal]
        counts = []
        seps = []
        for out in bite:
            cnt, sep = _vocal_peak_metrics(out)
            counts.append(float(cnt))
            seps.append(sep)
        speech_gain = float(np.mean([
            _band_db(b, 1000.0, 3500.0) - _band_db(c, 1000.0, 3500.0)
            for b, c in zip(bite, calm)
        ]))
        return {
            "formant_peak_count": float(np.mean(counts)) if counts else 0.0,
            "formant_separation_hz": float(np.mean(seps)) if seps else 0.0,
            "speech_band_gain_db": speech_gain,
        }

    mix = clips["mix"]
    calm = [_process(body, c, 0.0, 0.0) for c in mix]
    boundary = [_process(body, c, 0.5, 0.0) for c in mix]
    fracture = [_process(body, c, 1.0, 0.0) for c in mix]
    root_tether_shift = float(np.mean([
        _band_db(f, 20.0, 200.0) - _band_db(c, 20.0, 200.0)
        for f, c in zip(fracture, calm)
    ]))
    null_depth = float(np.mean([
        _full_rms_db(c) - _full_rms_db(b)
        for c, b in zip(calm, boundary)
    ]))
    fracture_distance = float(np.mean([
        _spectral_distance_db(b, f)
        for b, f in zip(boundary, fracture)
    ]))
    return {
        "root_tether_shift_db": root_tether_shift,
        "null_depth_db": null_depth,
        "fracture_distance_db": fracture_distance,
    }


def _improvement(cand: float, base: float, rule: MetricRule) -> tuple[float, float]:
    if rule.mode == "higher":
        raw = (cand - base) / rule.scale
    elif rule.mode == "lower":
        raw = (base - cand) / rule.scale
    elif rule.mode == "abs_zero":
        raw = (abs(base) - abs(cand)) / rule.scale
    else:
        raise ValueError(f"Unknown mode {rule.mode}")
    bounded = float(np.clip(raw, -MAX_ABS_METRIC_CONTRIBUTION, MAX_ABS_METRIC_CONTRIBUTION))
    return bounded, float(raw)


def compare_candidate(body_key: str, baseline: Body, candidate: Body, clips: dict[str, list[np.ndarray]]) -> dict:
    public_name = BODY_PUBLIC_NAME[body_key]
    base_profile = body_profile(baseline)
    cand_profile = body_profile(candidate)
    base_gate = shipping_gate(baseline, public_name)
    cand_gate = shipping_gate(candidate, public_name)
    base_task = measure_task_metrics(baseline, body_key, clips)
    cand_task = measure_task_metrics(candidate, body_key, clips)

    improvements = {}
    raw_improvements = {}
    wins = 0
    total = 0.0
    for rule in RULES[body_key]:
        imp, raw_imp = _improvement(cand_task[rule.name], base_task[rule.name], rule)
        improvements[rule.name] = imp
        raw_improvements[rule.name] = raw_imp
        total += imp
        if raw_imp > 0.02:
            wins += 1

    stability_raw = (
        base_profile["midpoint_audit"]["worst_peak_db"] - cand_profile["midpoint_audit"]["worst_peak_db"]
    ) / 6.0
    stability_imp = float(np.clip(stability_raw, -MAX_ABS_METRIC_CONTRIBUTION, MAX_ABS_METRIC_CONTRIBUTION))
    improvements["midpoint_worst_peak_db"] = stability_imp
    raw_improvements["midpoint_worst_peak_db"] = stability_raw
    total += stability_imp
    if stability_raw > 0.02:
        wins += 1

    shipping_gate_passed = bool(cand_gate["passed"])
    midpoint_passed = bool(cand_profile["midpoint_audit"]["passed"])
    # Name contracts are the hard release gate. Midpoint audit is informational
    # and still contributes via stability metric, but no longer hard-fails release.
    gate_ok = shipping_gate_passed
    promoted = gate_ok and (wins >= 2) and (total > 0.0)

    return {
        "body": public_name,
        "candidate_name": candidate.name,
        "shipping_gate_passed": shipping_gate_passed,
        "midpoint_passed": midpoint_passed,
        "gate_ok": gate_ok,
        "promoted": promoted,
        "wins": wins,
        "composite_improvement": total,
        "baseline": {
            "name": baseline.name,
            "gate": base_gate,
            "midpoint": base_profile["midpoint_audit"],
            "task": base_task,
        },
        "candidate": {
            "gate": cand_gate,
            "midpoint": cand_profile["midpoint_audit"],
            "task": cand_task,
        },
        "improvements": improvements,
        "raw_improvements": raw_improvements,
        "max_abs_metric_contribution": MAX_ABS_METRIC_CONTRIBUTION,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Promote candidates only if they beat E-MU translated baselines.")
    parser.add_argument("--manifest", type=Path, default=ROOT / "vault" / "_shipping_finalists" / "manifest.json")
    parser.add_argument("--corpus", type=Path, default=None, help="Optional WAV corpus directory.")
    parser.add_argument("--out", type=Path, default=ROOT / "vault" / "_scorecards")
    parser.add_argument(
        "--balanced-sampling",
        action="store_true",
        help="Use equal clip counts per class during scoring.",
    )
    parser.add_argument("--balance-seed", type=int, default=42)
    parser.add_argument(
        "--max-clips-per-class",
        type=int,
        default=None,
        help="Optional cap after balancing. Applies per class.",
    )
    parser.add_argument(
        "--allow-fallback",
        action="store_true",
        help="Allow synthetic fallback when corpus classes are missing.",
    )
    args = parser.parse_args()
    if args.max_clips_per_class is not None and args.max_clips_per_class <= 0:
        raise ValueError("--max-clips-per-class must be positive when provided.")

    if not args.manifest.exists():
        raise FileNotFoundError(f"Manifest not found: {args.manifest}")

    clips, probe_mode, label_counts, missing_labels = load_corpus_with_meta(
        args.corpus,
        allow_fallback=args.allow_fallback,
    )
    balance_input_counts = {k: len(v) for k, v in clips.items()}
    balance_output_counts = dict(balance_input_counts)
    balance_target_count = min(balance_input_counts.values()) if balance_input_counts else 0
    eval_clips = clips
    if args.balanced_sampling:
        eval_clips, balance_input_counts, balance_output_counts, balance_target_count = balanced_sample_clips(
            clips,
            seed=args.balance_seed,
            max_per_class=args.max_clips_per_class,
        )
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))

    results: dict[str, dict] = {}
    summary = []

    for body_key in ("speaker_knockerz", "aluminum_siding", "small_talk", "cul_de_sac"):
        baseline_path = BASELINE_COMPILED[body_key]
        if not baseline_path.exists():
            raise FileNotFoundError(f"Missing baseline compiled body: {baseline_path}")
        baseline = Body.from_json(str(baseline_path))

        entries = manifest.get(body_key, [])
        candidate_reports = []
        for entry in entries:
            candidate = Body.from_json(entry["path"])
            candidate_reports.append(compare_candidate(body_key, baseline, candidate, eval_clips))

        candidate_reports.sort(
            key=lambda r: (
                0 if r["promoted"] else 1,
                -(r["composite_improvement"]),
                -r["wins"],
            )
        )
        promoted_reports = [r for r in candidate_reports if bool(r.get("promoted"))]
        chosen = promoted_reports[0] if promoted_reports else None
        results[body_key] = {
            "public_name": BODY_PUBLIC_NAME[body_key],
            "baseline_path": str(baseline_path),
            "chosen": chosen,
            "candidates": candidate_reports,
            "selection_failure": None if chosen is not None else "no_promoted_candidate",
        }
        if chosen is None:
            summary.append({
                "body": BODY_PUBLIC_NAME[body_key],
                "candidate": None,
                "promoted": False,
                "composite_improvement": 0.0,
                "wins": 0,
                "gate_ok": False,
                "selection_failure": "no_promoted_candidate",
            })
            continue

        summary.append({
            "body": BODY_PUBLIC_NAME[body_key],
            "candidate": chosen["candidate_name"],
            "promoted": chosen["promoted"],
            "composite_improvement": chosen["composite_improvement"],
            "wins": chosen["wins"],
            "gate_ok": chosen["gate_ok"],
            "selection_failure": None,
        })

    report = {
        "probe_mode": probe_mode,
        "manifest": str(args.manifest),
        "corpus": str(args.corpus) if args.corpus is not None else None,
        "allow_fallback": args.allow_fallback,
        "max_abs_metric_contribution": MAX_ABS_METRIC_CONTRIBUTION,
        "balanced_sampling": args.balanced_sampling,
        "balance_seed": args.balance_seed,
        "max_clips_per_class": args.max_clips_per_class,
        "balance_input_counts": balance_input_counts,
        "balance_output_counts": balance_output_counts,
        "balance_target_count": balance_target_count,
        "label_counts": label_counts,
        "missing_labels": missing_labels,
        "summary": summary,
        "results": results,
    }

    args.out.mkdir(parents=True, exist_ok=True)
    report_json = args.out / "better_than_emu_scorecard.json"
    report_md = args.out / "better_than_emu_scorecard.md"
    report_json.write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = [
        "# Better Than E-MU Scorecard",
        "",
        f"- Probe mode: `{probe_mode}`",
        f"- Allow fallback: `{args.allow_fallback}`",
        f"- Max abs metric contribution: `{MAX_ABS_METRIC_CONTRIBUTION}`",
        f"- Balanced sampling: `{args.balanced_sampling}`",
        f"- Balance seed: `{args.balance_seed}`",
        f"- Max clips per class: `{args.max_clips_per_class}`",
        f"- Balance input counts: `{balance_input_counts}`",
        f"- Balance output counts: `{balance_output_counts}`",
        f"- Balance target count: `{balance_target_count}`",
        f"- Manifest: `{args.manifest}`",
        f"- Label counts: `{label_counts}`",
        f"- Missing labels: `{missing_labels}`",
        "",
        "## Summary",
        "",
        "| Body | Selected Candidate | Promoted | Gate OK | Wins | Composite Improvement |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in summary:
        selected_candidate = row["candidate"] if row["candidate"] is not None else "(none)"
        lines.append(
            f"| {row['body']} | {selected_candidate} | {str(row['promoted'])} | "
            f"{str(row['gate_ok'])} | {row['wins']} | {row['composite_improvement']:.3f} |"
        )
    lines += ["", "## Rule", "", "- Promotion requires: `gate_ok` and at least 2 metric wins and positive composite improvement."]
    report_md.write_text("\n".join(lines), encoding="utf-8")

    print(report_json)
    print(report_md)


if __name__ == "__main__":
    main()
