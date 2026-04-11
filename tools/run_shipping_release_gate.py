"""Fail-closed shipping release gate for emotional-vocal confidence.

Pipeline:
1) Build emotion-stress vocal corpus.
2) Create stratified design/holdout splits.
3) Run holdout promotion with balanced sampling.
4) Fail unless all shipping bodies are promoted and gate_ok on holdout.
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TOOLS = ROOT / "tools"
EXPECTED_BODIES = (
    "Speaker Knockerz",
    "Aluminum Siding",
    "Small Talk Ah-Ee",
    "Cul-De-Sac",
)


def _run(cmd: list[str]) -> None:
    print(f"$ {' '.join(cmd)}")
    subprocess.run(cmd, cwd=ROOT, check=True)


def _body_row_summary(row: dict) -> dict:
    return {
        "design_candidate": row.get("design_candidate"),
        "holdout_promoted": bool(row.get("holdout_promoted")),
        "holdout_gate_ok": bool(row.get("holdout_gate_ok")),
        "holdout_wins": int(row.get("holdout_wins", 0)),
        "holdout_composite_improvement": float(row.get("holdout_composite_improvement", 0.0)),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run fail-closed shipping release gate.")
    parser.add_argument(
        "--vocal-source",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "vocals_vocalset",
    )
    parser.add_argument(
        "--base-source",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "user_kits_real",
    )
    parser.add_argument(
        "--out-vocals",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "vocals_emotion_stress",
    )
    parser.add_argument(
        "--out-corpus",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "user_kits_plus_emotion_vocals",
    )
    parser.add_argument(
        "--design-corpus",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "design_split_emotion_strat",
    )
    parser.add_argument(
        "--holdout-corpus",
        type=Path,
        default=ROOT / "datasets" / "wav_corpus" / "holdout_split_emotion_strat",
    )
    parser.add_argument("--per-bucket", type=int, default=120)
    parser.add_argument("--max-per-class", type=int, default=70)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument(
        "--scorecard-out",
        type=Path,
        default=ROOT / "vault" / "_scorecards",
    )
    parser.add_argument(
        "--skip-corpus-build",
        action="store_true",
        help="Skip emotion corpus rebuild and use existing corpora/splits.",
    )
    args = parser.parse_args()

    python = sys.executable

    if not args.skip_corpus_build:
        _run(
            [
                python,
                str(TOOLS / "build_emotion_stress_corpus.py"),
                "--vocal-source",
                str(args.vocal_source),
                "--base-source",
                str(args.base_source),
                "--out-vocals",
                str(args.out_vocals),
                "--out-corpus",
                str(args.out_corpus),
                "--per-bucket",
                str(args.per_bucket),
                "--seed",
                str(args.seed),
            ]
        )
        _run(
            [
                python,
                str(TOOLS / "split_emotion_stress_corpus.py"),
                "--source",
                str(args.out_corpus),
                "--design",
                str(args.design_corpus),
                "--holdout",
                str(args.holdout_corpus),
                "--seed",
                str(args.seed),
                "--max-per-class",
                str(args.max_per_class),
                "--holdout-ratio",
                "0.2",
            ]
        )

    _run(
        [
            python,
            str(TOOLS / "promote_holdout.py"),
            "--design-corpus",
            str(args.design_corpus),
            "--holdout-corpus",
            str(args.holdout_corpus),
            "--out",
            str(args.scorecard_out),
            "--balanced-sampling",
            "--balance-seed",
            str(args.seed),
        ]
    )

    scorecard_path = args.scorecard_out / "better_than_emu_holdout_scorecard.json"
    if not scorecard_path.exists():
        raise FileNotFoundError(f"Expected scorecard not found: {scorecard_path}")
    scorecard = json.loads(scorecard_path.read_text(encoding="utf-8"))

    summary_rows = scorecard.get("summary", [])
    by_body = {row.get("body"): row for row in summary_rows}
    failures: list[str] = []
    selected: dict[str, dict] = {}

    for body in EXPECTED_BODIES:
        row = by_body.get(body)
        if row is None:
            failures.append(f"missing summary row for body '{body}'")
            continue
        selected[body] = _body_row_summary(row)
        if not row.get("holdout_gate_ok", False):
            failures.append(f"{body}: holdout_gate_ok=False")
        if not row.get("holdout_promoted", False):
            failures.append(f"{body}: holdout_promoted=False")
        if float(row.get("holdout_composite_improvement", 0.0)) <= 0.0:
            failures.append(f"{body}: holdout_composite_improvement<=0")
        if int(row.get("holdout_wins", 0)) < 2:
            failures.append(f"{body}: holdout_wins<2")

    gate_report = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "passed": len(failures) == 0,
        "expected_bodies": list(EXPECTED_BODIES),
        "failures": failures,
        "scorecard_path": str(scorecard_path),
        "protocol": scorecard.get("protocol", {}),
        "selected": selected,
    }

    args.scorecard_out.mkdir(parents=True, exist_ok=True)
    gate_json = args.scorecard_out / "shipping_release_gate.json"
    gate_md = args.scorecard_out / "shipping_release_gate.md"
    gate_json.write_text(json.dumps(gate_report, indent=2), encoding="utf-8")

    lines = [
        "# Shipping Release Gate",
        "",
        f"- Passed: `{gate_report['passed']}`",
        f"- Scorecard: `{scorecard_path}`",
        "",
        "## Bodies",
        "",
        "| Body | Candidate | Promoted | Gate OK | Wins | Composite |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for body in EXPECTED_BODIES:
        row = selected.get(body)
        if row is None:
            lines.append(f"| {body} | (missing) | False | False | 0 | 0.000 |")
            continue
        lines.append(
            f"| {body} | {row['design_candidate']} | {row['holdout_promoted']} | "
            f"{row['holdout_gate_ok']} | {row['holdout_wins']} | "
            f"{row['holdout_composite_improvement']:.3f} |"
        )

    if failures:
        lines += ["", "## Failures", ""]
        for item in failures:
            lines.append(f"- {item}")
    gate_md.write_text("\n".join(lines), encoding="utf-8")

    print(gate_json)
    print(gate_md)
    if failures:
        raise SystemExit("shipping release gate failed")


if __name__ == "__main__":
    main()
