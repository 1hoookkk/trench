#!/usr/bin/env python3
"""Plan a Bark frame map into heritage MorphDesigner authoring JSON.

Input:  trench.ui.frame_map.v1
Output: trench.heritage.morph_designer.v1, accepted by tools/compile_emu_designer.py

This is authoring/tooling only. It does not emit runtime coefficients.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any

try:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


AUTHORING_SR = 39062.5
NYQUIST_LIMIT_HZ = AUTHORING_SR * 0.5
VALID_RECIPES = {"talking_formant", "bass_pressure", "brittle_air", "peak_shelf"}


def bark_to_hz(bark: float) -> float:
    """Traunmueller inverse, matching tools/bark.py."""
    return 1960.0 * (bark + 0.53) / (26.28 - bark)


def hz_to_bark(hz: float) -> float:
    """Traunmueller forward, matching tools/bark.py."""
    return 26.81 * hz / (1960.0 + hz) - 0.53


class PlanError(Exception):
    pass


def _num(value: object, field: str) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError) as exc:
        raise PlanError(f"{field} must be numeric") from exc
    if not math.isfinite(out):
        raise PlanError(f"{field} must be finite")
    return out


def _freq(frame: dict[str, Any], landmark: str, frame_name: str) -> float:
    if landmark not in frame:
        raise PlanError(f"frames.{frame_name}.{landmark} missing")
    value = frame[landmark]
    if not isinstance(value, dict):
        raise PlanError(f"frames.{frame_name}.{landmark} must be an object")

    if "hz" in value:
        hz = _num(value["hz"], f"frames.{frame_name}.{landmark}.hz")
    elif "bark" in value:
        bark = _num(value["bark"], f"frames.{frame_name}.{landmark}.bark")
        hz = bark_to_hz(bark)
    else:
        raise PlanError(f"frames.{frame_name}.{landmark} needs hz or bark")

    if hz <= 0.0 or hz >= NYQUIST_LIMIT_HZ:
        raise PlanError(
            f"frames.{frame_name}.{landmark} frequency {hz:.3f} out of range"
        )
    return hz


def _gain(frame: dict[str, Any], landmark: str) -> float:
    value = frame.get(landmark)
    if isinstance(value, dict) and "gain_db" in value:
        return _num(value["gain_db"], f"{landmark}.gain_db")
    return 0.0


# Synthetic landmarks are computed from the frame map at planner time. They
# never appear in the frame_map fixture; the planner derives Hz from sibling
# landmarks on the same frame so the audit point stays in heritage Hz space.
#
# nasal_anti — the nasal antiformant sits between F1 and F2. The textbook
# location is the log-midpoint of the two formants (geometric mean in Hz,
# arithmetic mean in Bark), which gives a stable carve target whether the
# vowel is closed (F1≈342, F2≈2322) or open (F1≈700, F2≈1220).
SYNTHETIC_LANDMARKS: dict[str, dict[str, Any]] = {
    "nasal_anti": {
        "sources": ("throat_f1", "mouth_f2"),
        "method": "log_mid",
        "default_gain_db": -3.0,
    },
}


def _synthetic_present(
    frame_a: dict[str, Any], frame_b: dict[str, Any], name: str
) -> bool:
    spec = SYNTHETIC_LANDMARKS.get(name)
    if spec is None:
        return False
    return all(src in frame_a and src in frame_b for src in spec["sources"])


def _synthetic_freq(
    frame: dict[str, Any], name: str, frame_name: str
) -> float:
    spec = SYNTHETIC_LANDMARKS[name]
    sources = spec["sources"]
    method = spec["method"]
    freqs = [_freq(frame, src, frame_name) for src in sources]
    if method == "log_mid":
        product = 1.0
        for f in freqs:
            product *= f
        return product ** (1.0 / len(freqs))
    raise PlanError(f"unknown synthetic method {method!r} for {name}")


def _synthetic_gain(name: str) -> float:
    return float(SYNTHETIC_LANDMARKS[name]["default_gain_db"])


def _first_present(frame_a: dict[str, Any], frame_b: dict[str, Any], names: tuple[str, ...]) -> str:
    for name in names:
        if name in frame_a and name in frame_b:
            return name
        if _synthetic_present(frame_a, frame_b, name):
            return name
    raise PlanError(f"missing landmark pair: {'|'.join(names)}")


def load_frame_map(path: Path) -> dict[str, Any]:
    data = json.loads(path.read_text(encoding="utf-8"))
    if data.get("schema") != "trench.ui.frame_map.v1":
        raise PlanError(f"expected schema trench.ui.frame_map.v1, got {data.get('schema')}")
    frames = data.get("frames")
    if not isinstance(frames, dict) or not isinstance(frames.get("A"), dict) or not isinstance(frames.get("B"), dict):
        raise PlanError("missing frames.A or frames.B")
    return data


def _stage(
    shape: str,
    landmark_options: tuple[str, ...],
    frame_a: dict[str, Any],
    frame_b: dict[str, Any],
) -> dict[str, Any]:
    landmark = _first_present(frame_a, frame_b, landmark_options)
    if landmark in SYNTHETIC_LANDMARKS and landmark not in frame_a:
        fa = _synthetic_freq(frame_a, landmark, "A")
        fb = _synthetic_freq(frame_b, landmark, "B")
        gain_a = _synthetic_gain(landmark)
        gain_b = _synthetic_gain(landmark)
    else:
        fa = _freq(frame_a, landmark, "A")
        fb = _freq(frame_b, landmark, "B")
        gain_a = _gain(frame_a, landmark)
        gain_b = _gain(frame_b, landmark)
    return {
        "shape": shape,
        "lo": {"freq_hz": round(fa, 1), "gain_db": round(gain_a, 2)},
        "hi": {"freq_hz": round(fb, 1), "gain_db": round(gain_b, 2)},
        "planner_landmark": landmark,
    }


MORPH_RECIPES: dict[str, dict[str, Any]] = {
    "talking_formant": {
        "note": (
            "Heritage VOW skeleton: LP boundary / F1 throat / F2 mouth / "
            "nasal antiresonance / bite / HP cleanup. The S3 nasal cut is "
            "synthesised at log-mid(F1, F2) — the textbook antiformant zone "
            "between the throat and mouth peaks — so the carve never collides "
            "with bite or presence."
        ),
        "stages": (
            ("LP", ("air", "sibilance", "presence")),  # S0 — LP boundary / darkness
            ("EQ", ("throat_f1",)),                    # S1 — F1 throat formant
            ("EQ", ("mouth_f2",)),                     # S2 — F2 mouth formant
            ("ZERO", ("nasal_anti",)),                 # S3 — synthesised antiformant
            ("EQ", ("bite",)),                         # S4 — bite / edge injury
            ("HP", ("floor", "chest")),                # S5 — HP cleanup
        ),
    },
    "bass_pressure": {
        "note": "Heritage bass-pressure skeleton: LP + 4xEQ + HP.",
        "stages": (
            ("LP", ("air", "sibilance", "presence")),
            ("EQ", ("chest",)),
            ("EQ", ("choke", "throat_f1")),
            ("EQ", ("rip", "mouth_f2")),
            ("EQ", ("cry", "bite")),
            ("HP", ("floor",)),
        ),
    },
    "brittle_air": {
        "note": "Heritage brittle-air skeleton: LP + 4xEQ + HP.",
        "stages": (
            ("HP", ("floor",)),
            ("EQ", ("void", "chest")),
            ("EQ", ("bite",)),
            ("EQ", ("presence",)),
            ("EQ", ("sibilance", "air")),
            ("LP", ("air",)),
        ),
    },
}


def plan_morph_designer(frame_map: dict[str, Any], recipe: str) -> dict[str, Any]:
    frames = frame_map["frames"]
    frame_a = frames["A"]
    frame_b = frames["B"]
    spec = MORPH_RECIPES[recipe]
    boost = _num(frame_map.get("boost", 0.004), "boost")
    return {
        "schema": "trench.heritage.morph_designer.v1",
        "filter_class": "MorphDesigner",
        "name": frame_map.get("name", "unnamed"),
        "category": frame_map.get("category", "PROG"),
        "source_affinity": frame_map.get("source_affinity", []),
        "planner": {
            "input_schema": "trench.ui.frame_map.v1",
            "recipe": recipe,
            "recipe_note": spec["note"],
            "landmark_source": "authoring/ui/bark_landmarks.json",
            "authoring_sample_rate_hz": AUTHORING_SR,
        },
        "boost": boost,
        "stages": [
            _stage(shape, landmarks, frame_a, frame_b)
            for shape, landmarks in spec["stages"]
        ],
    }


def plan_peak_shelf(frame_map: dict[str, Any]) -> dict[str, Any]:
    frames = frame_map["frames"]

    def one(frame: dict[str, Any], frame_name: str) -> dict[str, Any]:
        landmark = _first_present(frame, frame, ("floor", "mouth_f2", "throat_f1"))
        freq = _freq(frame, landmark, frame_name)
        shelf = 0.0
        peak = 0.0
        if isinstance(frame.get("shelf"), dict):
            shelf = _num(frame["shelf"].get("value", 0.0), f"frames.{frame_name}.shelf.value")
        if isinstance(frame.get("peak_db"), dict):
            peak = _num(frame["peak_db"].get("value", 0.0), f"frames.{frame_name}.peak_db.value")
        return {"freq_hz": round(freq, 1), "shelf": shelf, "peak_db": peak}

    return {
        "schema": "trench.heritage.peak_shelf.v1",
        "filter_class": "PeakShelfMorph",
        "name": frame_map.get("name", "unnamed"),
        "category": frame_map.get("category", "PROG"),
        "source_affinity": frame_map.get("source_affinity", []),
        "planner": {
            "input_schema": "trench.ui.frame_map.v1",
            "recipe": "peak_shelf",
            "recipe_note": "Heritage Peak/Shelf Morph frame pair.",
        },
        "frames": {"A": one(frames["A"], "A"), "B": one(frames["B"], "B")},
    }


def validate_no_compiled(output: dict[str, Any]) -> list[str]:
    errors: list[str] = []

    def check(obj: Any, path: str = "") -> None:
        if isinstance(obj, dict):
            for key, value in obj.items():
                lower = key.lower()
                if lower in {
                    "m0_q0",
                    "m0_q100",
                    "m100_q0",
                    "m100_q100",
                    "coeffs",
                    "coefficient",
                    "compiled",
                    "keyframes",
                    "corners",
                }:
                    errors.append(f"{path}.{key}" if path else key)
                check(value, f"{path}.{key}" if path else key)
        elif isinstance(obj, list):
            for idx, value in enumerate(obj):
                check(value, f"{path}[{idx}]")

    check(output)
    return errors


def generate_plot(output: dict[str, Any], out_path: Path) -> None:
    if output.get("schema") != "trench.heritage.morph_designer.v1":
        raise PlanError("--plot only supports morph_designer output")
    stages = output.get("stages", [])
    labels = [f"S{i + 1}:{s.get('shape', '?')}" for i, s in enumerate(stages)]
    a_hz = [float(s["lo"]["freq_hz"]) for s in stages]
    b_hz = [float(s["hi"]["freq_hz"]) for s in stages]
    x = range(len(stages))

    fig, (ax_hz, ax_bark) = plt.subplots(2, 1, figsize=(9.5, 7.0), dpi=140)
    fig.suptitle(output.get("name", "unnamed"))

    ax_hz.plot(x, a_hz, marker="o", label="A / lo")
    ax_hz.plot(x, b_hz, marker="o", label="B / hi")
    ax_hz.set_yscale("log")
    ax_hz.set_ylabel("Hz")
    ax_hz.grid(True, which="both", alpha=0.3)
    ax_hz.legend()

    ax_bark.plot(x, [hz_to_bark(v) for v in a_hz], marker="o", label="A / lo")
    ax_bark.plot(x, [hz_to_bark(v) for v in b_hz], marker="o", label="B / hi")
    ax_bark.set_ylabel("Bark")
    ax_bark.grid(True, alpha=0.3)
    ax_bark.legend()

    ax_bark.set_xticks(list(x))
    ax_bark.set_xticklabels(labels, rotation=30, ha="right")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_path)
    plt.close(fig)
    print(out_path)


def print_summary(output: dict[str, Any]) -> None:
    print(f"schema={output.get('schema')}")
    print(f"recipe={output.get('planner', {}).get('recipe')}")
    if output.get("schema") == "trench.heritage.morph_designer.v1":
        for idx, stage in enumerate(output.get("stages", []), start=1):
            print(
                f"S{idx} {stage.get('shape')} {stage.get('planner_landmark')} "
                f"A={stage['lo']['freq_hz']}Hz B={stage['hi']['freq_hz']}Hz"
            )


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Stage planner: Bark Frame Map -> heritage authoring JSON"
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--recipe", required=True, choices=sorted(VALID_RECIPES))
    parser.add_argument("--out", required=True, type=Path)
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()

    try:
        frame_map = load_frame_map(args.input)
        if args.recipe == "peak_shelf":
            output = plan_peak_shelf(frame_map)
        else:
            output = plan_morph_designer(frame_map, args.recipe)
        errors = validate_no_compiled(output)
        if errors:
            raise PlanError(f"output contains compiled/runtime fields: {errors}")
    except (OSError, json.JSONDecodeError, PlanError) as exc:
        print(f"plan error: {exc}", file=sys.stderr)
        return 1

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(output, indent=2) + "\n", encoding="utf-8")
    print_summary(output)

    if args.plot:
        if not HAS_MATPLOTLIB:
            print("plot error: matplotlib not available", file=sys.stderr)
            return 1
        try:
            generate_plot(output, args.plot)
        except PlanError as exc:
            print(f"plot error: {exc}", file=sys.stderr)
            return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
