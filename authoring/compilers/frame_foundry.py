#!/usr/bin/env python3
"""Frame Foundry v0 candidate driver.

Reads a base MorphDesigner JSON + a mutation spec, produces N variants,
compiles, scores, rejects peaks > +45 dB, plots survivors (START | END),
and tiles a contact sheet.

Mutation spec shape:
{
  "candidates": [
    {
      "name": "A1_nasal_lift",
      "ops": [
        {"stage": 2, "frame": "hi", "freq_hz": 2400},
        {"stage": 3, "frame": "hi", "gain_db": 4.0}
      ]
    }
  ]
}

Only lo/hi.freq_hz and lo/hi.gain_db are mutable. Shape, stage count,
filter_class, and boost are frozen (they belong to the compiler surface).

Output layout under --out-dir:
  designers/ raw/ compiled/ scores/ plots/ rejects/
  contact_sheet.png   run_summary.json
"""
from __future__ import annotations

import argparse
import copy
import json
import math
import sys
from pathlib import Path

import numpy as np
import scipy.signal
import wave
import matplotlib
matplotlib.use("Agg")
import matplotlib.image as mpimg
import matplotlib.pyplot as plt

THIS = Path(__file__).resolve().parent
if str(THIS) not in sys.path:
    sys.path.insert(0, str(THIS))

import compile_emu_designer
import compile_raw
import plot_magnitude
import bark

HOST_SR = 44100.0
CORNERS = ("M0_Q0", "M0_Q100", "M100_Q0", "M100_Q100")
PASSTHROUGH_C = (1.0, 0.0, 0.0, 0.0, 0.0)
PEAK_REJECT_DB = 45.0
AUDITION_SECONDS = 4.0
AUDITION_SR = 44100
AUDITION_SEED = 42

BG = "#0b1020"
FG = "#cbd5f5"
EDGE = "#1e293b"


def _coefs(s: dict) -> tuple:
    return (float(s["c0"]), float(s["c1"]), float(s["c2"]),
            float(s["c3"]), float(s["c4"]))


def _cascade_db(stages: list[dict], boost: float, freqs: np.ndarray) -> np.ndarray:
    w = 2.0 * math.pi * freqs / HOST_SR
    e1, e2 = np.exp(-1j * w), np.exp(-2j * w)
    mag = np.ones_like(freqs)
    for s in stages:
        c = _coefs(s)
        if c == PASSTHROUGH_C:
            continue
        mag *= np.abs((c[0] + c[1]*e1 + c[2]*e2) / (1.0 + c[3]*e1 + c[4]*e2))
    return 20.0 * np.log10(np.maximum(mag * max(boost, 1e-9), 1e-12))


def _band_mean(freqs: np.ndarray, db: np.ndarray, lo: float, hi: float) -> float:
    mask = (freqs >= lo) & (freqs <= hi)
    return float(np.mean(db[mask])) if mask.any() else float("nan")


def _find_peaks(freqs: np.ndarray, db: np.ndarray, top_n: int = 6,
                min_db: float = 0.0) -> list[dict]:
    """Local-maxima peak finder. Reports freq_hz + bark + db per peak."""
    peaks: list[dict] = []
    for i in range(1, len(db) - 1):
        if db[i] > db[i - 1] and db[i] > db[i + 1] and db[i] > min_db:
            peaks.append({
                "freq_hz": float(freqs[i]),
                "bark": float(bark.hz_to_bark(freqs[i])),
                "db": float(db[i]),
            })
    peaks.sort(key=lambda p: p["db"], reverse=True)
    return peaks[:top_n]


def score_cartridge(cart: dict) -> dict:
    freqs = np.logspace(math.log10(20.0), math.log10(20000.0), 2048)
    peaks, sub, bright, active, peak_list = {}, {}, {}, {}, {}
    for kf in cart.get("keyframes", []):
        label = kf.get("label")
        if label not in CORNERS:
            continue
        stages = kf.get("stages", [])
        db = _cascade_db(stages, float(kf.get("boost", 1.0)), freqs)
        peaks[label] = float(np.max(db))
        sub[label] = _band_mean(freqs, db, 20.0, 60.0)  # Hz invariant: sub floor
        bright[label] = (_band_mean(freqs, db, 2000.0, 8000.0)
                         - _band_mean(freqs, db, 200.0, 800.0))
        active[label] = sum(1 for s in stages if _coefs(s) != PASSTHROUGH_C)
        peak_list[label] = _find_peaks(freqs, db)
    peak_max = max(peaks.values()) if peaks else float("nan")
    sb, eb = bright.get("M0_Q0", float("nan")), bright.get("M100_Q0", float("nan"))

    # Bark-space morph-delta per corner pair: how far (in Bark) do the
    # dominant peaks move from M0 to M100? Frequency target error is a
    # Bark distance; Hz invariants (sub_band_20_60_db) are kept above.
    def _top_bark(label):
        lst = peak_list.get(label, [])
        return lst[0]["bark"] if lst else float("nan")
    peak_bark_delta = {
        "top_peak_bark_m0_to_m100": float(_top_bark("M100_Q0") - _top_bark("M0_Q0"))
            if (peak_list.get("M0_Q0") and peak_list.get("M100_Q0")) else float("nan"),
    }

    return {
        "name": cart.get("name", "unnamed"),
        "peak_db_max": peak_max,
        "peak_db_per_corner": peaks,
        "sub_band_20_60_db": sub,
        "brightness_db": bright,
        "end_brighter_than_start": math.isfinite(sb) and math.isfinite(eb) and eb > sb,
        "active_stage_count": active,
        "peaks": peak_list,
        "peak_bark_delta": peak_bark_delta,
    }


def apply_ops(base: dict, ops: list[dict], name: str) -> dict:
    variant = copy.deepcopy(base)
    stages = variant.get("stages", [])
    for i, op in enumerate(ops):
        path = f"{name}.ops[{i}]"
        idx, frame = op.get("stage"), op.get("frame")
        if not isinstance(idx, int) or idx < 0 or idx >= len(stages):
            raise ValueError(f"{path}: stage index out of range")
        if frame not in ("lo", "hi"):
            raise ValueError(f"{path}: frame must be 'lo' or 'hi'")
        keys = [k for k in ("freq_hz", "gain_db") if k in op]
        if len(keys) != 1:
            raise ValueError(f"{path}: set exactly one of freq_hz or gain_db")
        stages[idx].setdefault(frame, {})[keys[0]] = float(op[keys[0]])
    return variant


def _pink_noise(n: int, rng: np.random.Generator) -> np.ndarray:
    """Paul Kellet's IIR pink-noise approximation. Deterministic via seeded rng."""
    white = rng.standard_normal(n)
    b = [0.049922035, -0.095993537, 0.050612699, -0.004408786]
    a = [1.0, -2.494956002, 2.017265875, -0.522189400]
    pink = scipy.signal.lfilter(b, a, white)
    peak = float(np.max(np.abs(pink)))
    return pink / peak * 0.6 if peak > 0 else pink


def _cascade_render(signal: np.ndarray, stages: list[dict], boost: float) -> np.ndarray:
    y = signal.astype(np.float64)
    for s in stages:
        c = _coefs(s)
        if c == PASSTHROUGH_C:
            continue
        y = scipy.signal.lfilter([c[0], c[1], c[2]], [1.0, c[3], c[4]], y)
    return y * max(boost, 1e-9)


def _find_corner(cart: dict, label: str) -> dict:
    for kf in cart.get("keyframes", []):
        if kf.get("label") == label:
            return kf
    raise ValueError(f"corner {label} not found")


def render_audition(cart: dict, out_path: Path,
                    duration_s: float = AUDITION_SECONDS,
                    sr: int = AUDITION_SR) -> Path:
    """Pink-noise morph audition: M0 -> M100 linear crossfade.

    Signal-side crossfade (approximation of runtime per-sample coefficient
    interpolation). Deterministic via seeded rng. Mono 16-bit PCM.
    """
    n = int(duration_s * sr)
    rng = np.random.default_rng(AUDITION_SEED)
    src = _pink_noise(n, rng)

    m0 = _find_corner(cart, "M0_Q0")
    m100 = _find_corner(cart, "M100_Q0")
    y0 = _cascade_render(src, m0["stages"], float(m0.get("boost", 1.0)))
    y1 = _cascade_render(src, m100["stages"], float(m100.get("boost", 1.0)))

    alpha = np.linspace(0.0, 1.0, n)
    y = y0 * (1.0 - alpha) + y1 * alpha

    peak = float(np.max(np.abs(y)))
    if peak > 1e-9:
        y = y * (0.92 / peak)

    pcm = np.clip(y * 32767.0, -32768, 32767).astype(np.int16)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with wave.open(str(out_path), "wb") as w:
        w.setnchannels(1)
        w.setsampwidth(2)
        w.setframerate(sr)
        w.writeframes(pcm.tobytes())
    return out_path


def _caption(score: dict) -> str:
    peak = score["peak_db_max"]
    bs = score["brightness_db"].get("M0_Q0", float("nan"))
    be = score["brightness_db"].get("M100_Q0", float("nan"))
    arrow = "bright+" if score["end_brighter_than_start"] else "bright-"
    return (f"{score['name']}\n"
            f"peak {peak:+.1f} dB   bright {bs:+.1f} -> {be:+.1f}   {arrow}")


def _build_sheet(run_dir: Path, out_path: Path) -> Path:
    names = sorted(p.stem for p in (run_dir / "plots").glob("*.png"))
    if not names:
        raise SystemExit(f"no plots in {run_dir / 'plots'}")
    cols = min(4, len(names))
    rows = math.ceil(len(names) / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(cols * 4.0, rows * 3.2),
                             dpi=130, facecolor=BG)
    axes_arr = np.atleast_1d(axes).flatten()
    for i, name in enumerate(names):
        ax = axes_arr[i]
        ax.imshow(mpimg.imread(run_dir / "plots" / f"{name}.png"))
        ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_color(EDGE)
        score_path = run_dir / "scores" / f"{name}.json"
        caption = (_caption(json.loads(score_path.read_text(encoding="utf-8")))
                   if score_path.exists() else name)
        ax.set_title(caption, color=FG, fontsize=8, loc="left", pad=4)
    for ax in axes_arr[len(names):]:
        ax.set_visible(False)
    fig.tight_layout()
    fig.savefig(out_path, facecolor=BG)
    plt.close(fig)
    return out_path


def run(base_path: Path, mutations_path: Path, out_dir: Path,
        plot_xaxis: str = "hz") -> dict:
    for sub in ("designers", "raw", "compiled", "plots", "scores", "rejects", "audition"):
        (out_dir / sub).mkdir(parents=True, exist_ok=True)

    base = json.loads(base_path.read_text(encoding="utf-8"))
    muts = json.loads(mutations_path.read_text(encoding="utf-8"))

    accepted: list[str] = []
    rejected: list[dict] = []

    for cand in muts.get("candidates", []):
        name = cand["name"]
        variant = apply_ops(base, cand.get("ops", []), name)
        variant["name"] = name
        (out_dir / "designers" / f"{name}.json").write_text(
            json.dumps(variant, indent=2) + "\n", encoding="utf-8")

        try:
            raw = compile_emu_designer.compile_designer(variant)
        except compile_emu_designer.CompileError as exc:
            rejected.append({"name": name, "stage": "designer", "reason": str(exc)})
            (out_dir / "rejects" / f"{name}.reason.txt").write_text(
                f"designer compile: {exc}\n", encoding="utf-8")
            continue
        (out_dir / "raw" / f"{name}.raw.json").write_text(
            json.dumps(raw, indent=2) + "\n", encoding="utf-8")

        try:
            cart = compile_raw.compile_raw(raw)
        except compile_raw.CompileError as exc:
            rejected.append({"name": name, "stage": "raw", "reason": str(exc)})
            (out_dir / "rejects" / f"{name}.reason.txt").write_text(
                f"raw compile: {exc}\n", encoding="utf-8")
            continue
        (out_dir / "compiled" / f"{name}.compiled.json").write_text(
            json.dumps(cart, indent=2) + "\n", encoding="utf-8")

        result = score_cartridge(cart)
        (out_dir / "scores" / f"{name}.json").write_text(
            json.dumps(result, indent=2) + "\n", encoding="utf-8")

        peak = result["peak_db_max"]
        if peak > PEAK_REJECT_DB:
            reason = f"peak_db_max {peak:+.1f} > {PEAK_REJECT_DB:+.1f}"
            rejected.append({"name": name, "stage": "score", "reason": reason})
            (out_dir / "rejects" / f"{name}.reason.txt").write_text(reason + "\n",
                                                                    encoding="utf-8")
            continue

        plot_magnitude.plot(cart, out_dir / "plots" / f"{name}.png",
                            selected=["M0_Q0", "M100_Q0"],
                            xaxis=plot_xaxis)
        render_audition(cart, out_dir / "audition" / f"{name}.wav")
        accepted.append(name)

    sheet_path = _build_sheet(out_dir, out_dir / "contact_sheet.png") if accepted else None
    summary = {
        "base": str(base_path),
        "mutations": str(mutations_path),
        "out_dir": str(out_dir),
        "thresholds": {"peak_reject_db": PEAK_REJECT_DB},
        "accepted": accepted,
        "rejected": rejected,
        "contact_sheet": str(sheet_path) if sheet_path else None,
    }
    (out_dir / "run_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    return summary


def main(argv: list[str]) -> int:
    ap = argparse.ArgumentParser(description="Frame Foundry candidate driver.")
    ap.add_argument("base", type=Path)
    ap.add_argument("--mutations", type=Path, required=True)
    ap.add_argument("--out-dir", type=Path, required=True)
    ap.add_argument("--plot-xaxis", choices=("hz", "bark"), default="hz",
                    help="X-axis for per-candidate plots (default: hz).")
    args = ap.parse_args(argv)

    summary = run(args.base, args.mutations, args.out_dir,
                  plot_xaxis=args.plot_xaxis)
    print(json.dumps({
        "accepted": len(summary["accepted"]),
        "rejected": len(summary["rejected"]),
        "contact_sheet": summary["contact_sheet"],
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
