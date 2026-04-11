"""Holdout promotion protocol for "better than E-MU".

Protocol:
1. Rank candidates on DESIGN corpus only.
2. Freeze the top design candidate per body.
3. Evaluate frozen winners on HOLDOUT corpus.
4. Promote only from holdout result.

This prevents tuning and judging on the same data.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
import sys
sys.path.insert(0, str(ROOT))

import tools.promote_better_than_emu as scorer
from pyruntime.body import Body

BODY_KEYS = ("speaker_knockerz", "aluminum_siding", "small_talk", "cul_de_sac")


def _validate_corpus_exists(corpus_dir: Path) -> None:
    if not corpus_dir.exists():
        raise FileNotFoundError(f"Corpus directory not found: {corpus_dir}")


def _design_rank_key(report: dict) -> tuple[int, float, int]:
    return (
        0 if report["gate_ok"] else 1,
        -report["composite_improvement"],
        -report["wins"],
    )


def main() -> None:
    parser = argparse.ArgumentParser(description="Design/holdout promotion protocol.")
    parser.add_argument("--manifest", type=Path, default=ROOT / "vault" / "_shipping_finalists" / "manifest.json")
    parser.add_argument("--design-corpus", type=Path, required=True)
    parser.add_argument("--holdout-corpus", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=ROOT / "vault" / "_scorecards")
    parser.add_argument(
        "--balanced-sampling",
        action="store_true",
        help="Use equal clip counts per class during both design and holdout evaluation.",
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

    _validate_corpus_exists(args.design_corpus)
    _validate_corpus_exists(args.holdout_corpus)
    design_clips, design_mode, design_counts, design_missing = scorer.load_corpus_with_meta(
        args.design_corpus,
        allow_fallback=args.allow_fallback,
    )
    holdout_clips, holdout_mode, holdout_counts, holdout_missing = scorer.load_corpus_with_meta(
        args.holdout_corpus,
        allow_fallback=args.allow_fallback,
    )
    design_balance_input = {k: len(v) for k, v in design_clips.items()}
    design_balance_output = dict(design_balance_input)
    design_balance_target = min(design_balance_input.values()) if design_balance_input else 0
    holdout_balance_input = {k: len(v) for k, v in holdout_clips.items()}
    holdout_balance_output = dict(holdout_balance_input)
    holdout_balance_target = min(holdout_balance_input.values()) if holdout_balance_input else 0
    if args.balanced_sampling:
        design_clips, design_balance_input, design_balance_output, design_balance_target = scorer.balanced_sample_clips(
            design_clips,
            seed=args.balance_seed,
            max_per_class=args.max_clips_per_class,
        )
        holdout_clips, holdout_balance_input, holdout_balance_output, holdout_balance_target = scorer.balanced_sample_clips(
            holdout_clips,
            seed=args.balance_seed + 1,
            max_per_class=args.max_clips_per_class,
        )

    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    results: dict[str, dict] = {}
    summary: list[dict] = []

    for body_key in BODY_KEYS:
        baseline_path = scorer.BASELINE_COMPILED[body_key]
        if not baseline_path.exists():
            raise FileNotFoundError(f"Missing baseline body: {baseline_path}")
        baseline = Body.from_json(str(baseline_path))

        entries = manifest.get(body_key, [])
        if not entries:
            results[body_key] = {
                "public_name": scorer.BODY_PUBLIC_NAME[body_key],
                "baseline_path": str(baseline_path),
                "design_ranked": [],
                "chosen_from_design": None,
                "holdout_eval": None,
                "holdout_oracle_top": None,
            }
            continue

        design_reports = []
        holdout_reports = []
        for entry in entries:
            candidate = Body.from_json(entry["path"])
            design_report = scorer.compare_candidate(body_key, baseline, candidate, design_clips)
            holdout_report = scorer.compare_candidate(body_key, baseline, candidate, holdout_clips)
            design_reports.append(design_report)
            holdout_reports.append(holdout_report)

        design_reports.sort(key=_design_rank_key)
        holdout_oracle = sorted(holdout_reports, key=_design_rank_key)[0]

        gate_passing_design = [r for r in design_reports if bool(r.get("gate_ok"))]
        chosen_design = gate_passing_design[0] if gate_passing_design else None

        if chosen_design is None:
            failure_reason = "no_gate_passing_design_candidate"
            results[body_key] = {
                "public_name": scorer.BODY_PUBLIC_NAME[body_key],
                "baseline_path": str(baseline_path),
                "design_ranked": design_reports,
                "chosen_from_design": None,
                "holdout_eval": None,
                "holdout_oracle_top": holdout_oracle,
                "selection_failure": failure_reason,
            }
            summary.append({
                "body": scorer.BODY_PUBLIC_NAME[body_key],
                "design_candidate": None,
                "design_gate_ok": False,
                "design_composite_improvement": 0.0,
                "holdout_promoted": False,
                "holdout_gate_ok": False,
                "holdout_wins": 0,
                "holdout_composite_improvement": 0.0,
                "holdout_oracle_candidate": holdout_oracle["candidate_name"],
                "holdout_oracle_promoted": holdout_oracle["promoted"],
                "selection_failure": failure_reason,
            })
            continue

        chosen_name = chosen_design["candidate_name"]
        chosen_holdout = next(r for r in holdout_reports if r["candidate_name"] == chosen_name)

        results[body_key] = {
            "public_name": scorer.BODY_PUBLIC_NAME[body_key],
            "baseline_path": str(baseline_path),
            "design_ranked": design_reports,
            "chosen_from_design": chosen_design,
            "holdout_eval": chosen_holdout,
            "holdout_oracle_top": holdout_oracle,
            "selection_failure": None,
        }
        summary.append({
            "body": scorer.BODY_PUBLIC_NAME[body_key],
            "design_candidate": chosen_name,
            "design_gate_ok": chosen_design["gate_ok"],
            "design_composite_improvement": chosen_design["composite_improvement"],
            "holdout_promoted": chosen_holdout["promoted"],
            "holdout_gate_ok": chosen_holdout["gate_ok"],
            "holdout_wins": chosen_holdout["wins"],
            "holdout_composite_improvement": chosen_holdout["composite_improvement"],
            "holdout_oracle_candidate": holdout_oracle["candidate_name"],
            "holdout_oracle_promoted": holdout_oracle["promoted"],
            "selection_failure": None,
        })

    report = {
        "protocol": {
            "manifest": str(args.manifest),
            "design_corpus": str(args.design_corpus),
            "holdout_corpus": str(args.holdout_corpus),
            "design_probe_mode": design_mode,
            "holdout_probe_mode": holdout_mode,
            "design_label_counts": design_counts,
            "holdout_label_counts": holdout_counts,
            "design_missing_labels": design_missing,
            "holdout_missing_labels": holdout_missing,
            "allow_fallback": args.allow_fallback,
            "max_abs_metric_contribution": scorer.MAX_ABS_METRIC_CONTRIBUTION,
            "balanced_sampling": args.balanced_sampling,
            "balance_seed": args.balance_seed,
            "max_clips_per_class": args.max_clips_per_class,
            "design_balance_input_counts": design_balance_input,
            "design_balance_output_counts": design_balance_output,
            "design_balance_target_count": design_balance_target,
            "holdout_balance_input_counts": holdout_balance_input,
            "holdout_balance_output_counts": holdout_balance_output,
            "holdout_balance_target_count": holdout_balance_target,
            "selection_rule": "Pick top candidate on design only; promote from holdout only.",
        },
        "summary": summary,
        "results": results,
    }

    args.out.mkdir(parents=True, exist_ok=True)
    report_json = args.out / "better_than_emu_holdout_scorecard.json"
    report_md = args.out / "better_than_emu_holdout_scorecard.md"
    report_json.write_text(json.dumps(report, indent=2), encoding="utf-8")

    lines = [
        "# Better Than E-MU Holdout Scorecard",
        "",
        f"- Design corpus: `{args.design_corpus}`",
        f"- Holdout corpus: `{args.holdout_corpus}`",
        f"- Design probe mode: `{design_mode}`",
        f"- Holdout probe mode: `{holdout_mode}`",
        f"- Allow fallback: `{args.allow_fallback}`",
        f"- Max abs metric contribution: `{scorer.MAX_ABS_METRIC_CONTRIBUTION}`",
        f"- Balanced sampling: `{args.balanced_sampling}`",
        f"- Balance seed: `{args.balance_seed}`",
        f"- Max clips per class: `{args.max_clips_per_class}`",
        f"- Design balance input counts: `{design_balance_input}`",
        f"- Design balance output counts: `{design_balance_output}`",
        f"- Design balance target count: `{design_balance_target}`",
        f"- Holdout balance input counts: `{holdout_balance_input}`",
        f"- Holdout balance output counts: `{holdout_balance_output}`",
        f"- Holdout balance target count: `{holdout_balance_target}`",
        f"- Design label counts: `{design_counts}`",
        f"- Holdout label counts: `{holdout_counts}`",
        f"- Design missing labels: `{design_missing}`",
        f"- Holdout missing labels: `{holdout_missing}`",
        "",
        "## Summary",
        "",
        "| Body | Design Winner | Holdout Promoted | Holdout Gate OK | Holdout Wins | Holdout Composite | Holdout Oracle |",
        "|---|---|---:|---:|---:|---:|---|",
    ]
    for row in summary:
        design_candidate = row["design_candidate"] if row["design_candidate"] is not None else "(none)"
        lines.append(
            f"| {row['body']} | {design_candidate} | {row['holdout_promoted']} | "
            f"{row['holdout_gate_ok']} | {row['holdout_wins']} | "
            f"{row['holdout_composite_improvement']:.3f} | {row['holdout_oracle_candidate']} |"
        )
    lines += [
        "",
        "## Rule",
        "",
        "- Candidate selection is based on design corpus only.",
        "- Final promotion decision is holdout-only.",
    ]
    report_md.write_text("\n".join(lines), encoding="utf-8")

    print(report_json)
    print(report_md)


if __name__ == "__main__":
    main()
