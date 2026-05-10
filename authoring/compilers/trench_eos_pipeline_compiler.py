#!/usr/bin/env python3
"""
TRENCH / E-MU EOS Pipeline Compiler Harness
==========================================

Internal reference-surface harness for the verified EOS-style runtime shape:

    1. Compiler      : 6-byte stage descriptors -> encoded uint16 quintuplets
    2. Interpolator  : bilinear interpolation in encoded c-domain
    3. Decoder       : uint16 quintuplets -> float coefficients c0..c4
    4. Render kernel : proprietary Rossum 6-stage serial biquad cascade

Boundary:
This file does not invent E-MU descriptor lowering, coefficient encoding,
decoder behavior, or the render equation. Those must come from verified reverse
engineering. Anything else would be fake precision.

What this does now:
- Defines the runtime data model.
- Preserves encoded uint16 coefficient space as the interpolation domain.
- Models Morph Designer as 2 logical authored corners duplicated into 4 runtime
  banks.
- Writes safe JSON with revision/hash/backup/atomic replace.
- Provides CLI and Tkinter wrappers for compiling/saving encoded fixtures.
- Exposes explicit replacement seams for the real lowerer/decoder/kernel.

Use this as a stable shell around the real lowering code. It is not a playable
compiled-v1 cartridge writer.
"""

from __future__ import annotations

import argparse
import datetime as _dt
import hashlib
import json
import os
import shutil
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, List, Mapping, Optional, Sequence, Tuple


COMPILER_NAME = "trench_eos_pipeline_compiler"
COMPILER_VERSION = "1.0.1"
SCHEMA = "trench.eos.encoded_surface.v1"

NATIVE_SAMPLE_RATE_HZ = 39_062.5
STAGE_COUNT = 6
COEFFS_PER_STAGE = 5
DESCRIPTOR_BYTES = 6
UINT16_MIN = 0
UINT16_MAX = 0xFFFF

DEFAULT_OUT_PATH = Path("cartridges") / "factory" / "generated" / "live_eos_surface.json"
DEFAULT_FIXTURE_PATH = Path("cartridges") / "factory" / "fixtures" / "encoded_fixture.json"

RUNTIME_CORNER_ORDER = ("m0_q0", "m1_q0", "m0_q1", "m1_q1")


def now_iso_utc() -> str:
    return _dt.datetime.now(_dt.timezone.utc).replace(microsecond=0).isoformat()


def stable_json_dumps(data: Any) -> str:
    return json.dumps(data, indent=2, sort_keys=True, ensure_ascii=False)


def short_hash(data: Mapping[str, Any]) -> str:
    return hashlib.sha256(stable_json_dumps(data).encode("utf-8")).hexdigest()[:16]


def project_root_from_script() -> Path:
    return Path(__file__).resolve().parents[2]


def resolve_path(path: Path) -> Path:
    if path.is_absolute():
        return path
    return project_root_from_script() / path


def clamp_int(value: int, lo: int, hi: int) -> int:
    return max(lo, min(hi, int(value)))


def validate_uint8(value: int, label: str) -> int:
    if not isinstance(value, int):
        raise TypeError(f"{label}: expected int uint8, got {type(value).__name__}")
    if not 0 <= value <= 255:
        raise ValueError(f"{label}: uint8 out of range: {value}")
    return value


def validate_uint16(value: int, label: str) -> int:
    if not isinstance(value, int):
        raise TypeError(f"{label}: expected int uint16, got {type(value).__name__}")
    if not UINT16_MIN <= value <= UINT16_MAX:
        raise ValueError(f"{label}: uint16 out of range: {value}")
    return value


def lerp_u16(a: int, b: int, t: float) -> int:
    """
    Runtime-like c-domain LERP for encoded coefficient words.

    Exact EOS rounding may use integer/fixed-point details. This uses nearest
    integer after scalar interpolation until the decompiled rounding path is
    wired in.
    """
    t = max(0.0, min(1.0, float(t)))
    return clamp_int(round(a + (b - a) * t), UINT16_MIN, UINT16_MAX)


def sqrt_compensation_level(host_sample_rate_hz: float) -> int:
    """
    Return the inferred EOS compensation level.

    Exact branch thresholds must be verified against the binary. This only
    classifies obvious integer-octave-rate cases.
    """
    ratio = float(host_sample_rate_hz) / NATIVE_SAMPLE_RATE_HZ
    if ratio < 1.5:
        return 0
    if ratio < 3.0:
        return 1
    return 2


@dataclass(frozen=True)
class StageDescriptor6:
    """One raw 6-byte stage descriptor from the EOS compiler input space."""

    raw: bytes

    def __post_init__(self) -> None:
        if len(self.raw) != DESCRIPTOR_BYTES:
            raise ValueError(
                f"StageDescriptor6 must be exactly {DESCRIPTOR_BYTES} bytes, got {len(self.raw)}"
            )

    @classmethod
    def from_hex(cls, hex_string: str) -> "StageDescriptor6":
        clean = hex_string.strip().replace(" ", "").replace("_", "")
        if len(clean) != DESCRIPTOR_BYTES * 2:
            raise ValueError(f"descriptor hex must be {DESCRIPTOR_BYTES * 2} hex chars")
        return cls(bytes.fromhex(clean))

    @classmethod
    def from_ints(cls, values: Sequence[int]) -> "StageDescriptor6":
        if len(values) != DESCRIPTOR_BYTES:
            raise ValueError(f"descriptor must contain {DESCRIPTOR_BYTES} byte values")
        raw = bytes(validate_uint8(v, f"descriptor_byte[{i}]") for i, v in enumerate(values))
        return cls(raw)

    def to_hex(self) -> str:
        return self.raw.hex()

    def to_byte_list(self) -> List[int]:
        return list(self.raw)


@dataclass(frozen=True)
class EncodedQuintuplet:
    """Five encoded uint16 coefficient words for one runtime stage."""

    words: Tuple[int, int, int, int, int]

    def __post_init__(self) -> None:
        if len(self.words) != COEFFS_PER_STAGE:
            raise ValueError(f"encoded quintuplet must have {COEFFS_PER_STAGE} words")
        for i, word in enumerate(self.words):
            validate_uint16(word, f"encoded_word[{i}]")

    @classmethod
    def from_list(cls, values: Sequence[int]) -> "EncodedQuintuplet":
        if len(values) != COEFFS_PER_STAGE:
            raise ValueError(f"encoded quintuplet must have {COEFFS_PER_STAGE} words")
        words = tuple(validate_uint16(int(v), f"encoded_word[{i}]") for i, v in enumerate(values))
        return cls(words)  # type: ignore[arg-type]

    @classmethod
    def identity_unknown(cls) -> "EncodedQuintuplet":
        raise NotImplementedError("No safe identity encoded quintuplet without verified EOS decoder/kernel")

    def to_list(self) -> List[int]:
        return list(self.words)

    def lerp(self, other: "EncodedQuintuplet", t: float) -> "EncodedQuintuplet":
        words = tuple(lerp_u16(a, b, t) for a, b in zip(self.words, other.words))
        return EncodedQuintuplet(words)  # type: ignore[arg-type]


@dataclass(frozen=True)
class EncodedStageBank:
    """Six serial stages, each represented by five encoded uint16 words."""

    stages: Tuple[
        EncodedQuintuplet,
        EncodedQuintuplet,
        EncodedQuintuplet,
        EncodedQuintuplet,
        EncodedQuintuplet,
        EncodedQuintuplet,
    ]

    def __post_init__(self) -> None:
        if len(self.stages) != STAGE_COUNT:
            raise ValueError(f"stage bank must contain exactly {STAGE_COUNT} stages")

    @classmethod
    def from_matrix(cls, matrix: Sequence[Sequence[int]]) -> "EncodedStageBank":
        if len(matrix) != STAGE_COUNT:
            raise ValueError(f"stage bank matrix must contain {STAGE_COUNT} rows")
        stages = tuple(EncodedQuintuplet.from_list(row) for row in matrix)
        return cls(stages)  # type: ignore[arg-type]

    def to_matrix(self) -> List[List[int]]:
        return [stage.to_list() for stage in self.stages]

    def lerp(self, other: "EncodedStageBank", t: float) -> "EncodedStageBank":
        stages = tuple(a.lerp(b, t) for a, b in zip(self.stages, other.stages))
        return EncodedStageBank(stages)  # type: ignore[arg-type]


@dataclass(frozen=True)
class RuntimeEncodedSurface:
    """
    Four runtime banks used by EOS interpolation.

    Morph Designer authored only two logical corners. In that case:
      m0_q0 == m0_q1 == logical m0
      m1_q0 == m1_q1 == logical m1
    """

    m0_q0: EncodedStageBank
    m1_q0: EncodedStageBank
    m0_q1: EncodedStageBank
    m1_q1: EncodedStageBank

    @classmethod
    def from_morph_designer_pair(cls, m0: EncodedStageBank, m1: EncodedStageBank) -> "RuntimeEncodedSurface":
        return cls(m0_q0=m0, m1_q0=m1, m0_q1=m0, m1_q1=m1)

    @classmethod
    def from_runtime_corners(cls, corners: Mapping[str, EncodedStageBank]) -> "RuntimeEncodedSurface":
        missing = [name for name in RUNTIME_CORNER_ORDER if name not in corners]
        if missing:
            raise ValueError(f"missing runtime corner(s): {missing}")
        return cls(
            m0_q0=corners["m0_q0"],
            m1_q0=corners["m1_q0"],
            m0_q1=corners["m0_q1"],
            m1_q1=corners["m1_q1"],
        )

    def to_json_corners(self) -> Dict[str, List[List[int]]]:
        return {
            "m0_q0": self.m0_q0.to_matrix(),
            "m1_q0": self.m1_q0.to_matrix(),
            "m0_q1": self.m0_q1.to_matrix(),
            "m1_q1": self.m1_q1.to_matrix(),
        }

    def interpolate_encoded(self, morph: float, q: float) -> EncodedStageBank:
        """Bilinear interpolation in encoded uint16 c-domain."""
        morph = max(0.0, min(1.0, float(morph)))
        q = max(0.0, min(1.0, float(q)))
        bottom = self.m0_q0.lerp(self.m1_q0, morph)
        top = self.m0_q1.lerp(self.m1_q1, morph)
        return bottom.lerp(top, q)


@dataclass(frozen=True)
class LogicalDescriptorBody:
    """Morph Designer style body: two logical corners, each with six descriptors."""

    name: str
    m0: Tuple[
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
    ]
    m1: Tuple[
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
        StageDescriptor6,
    ]

    @classmethod
    def from_json_dict(cls, data: Mapping[str, Any]) -> "LogicalDescriptorBody":
        name = str(data.get("name", "Untitled_EOS_Body"))
        logical = data.get("logical_descriptors") or data.get("descriptors")
        if not isinstance(logical, Mapping):
            raise ValueError("fixture must contain logical_descriptors or descriptors")

        def read_corner(corner: str) -> Tuple[StageDescriptor6, ...]:
            rows = logical.get(corner)
            if not isinstance(rows, list) or len(rows) != STAGE_COUNT:
                raise ValueError(f"logical corner {corner!r} must contain {STAGE_COUNT} descriptors")
            descs: List[StageDescriptor6] = []
            for i, row in enumerate(rows):
                if isinstance(row, str):
                    descs.append(StageDescriptor6.from_hex(row))
                elif isinstance(row, list):
                    descs.append(StageDescriptor6.from_ints(row))
                else:
                    raise TypeError(f"descriptor {corner}[{i}] must be hex string or list of six bytes")
            return tuple(descs)

        return cls(name=name, m0=read_corner("m0"), m1=read_corner("m1"))  # type: ignore[arg-type]

    def to_json_dict(self) -> Dict[str, Any]:
        return {
            "name": self.name,
            "logical_descriptors": {
                "m0": [d.to_hex() for d in self.m0],
                "m1": [d.to_hex() for d in self.m1],
            },
        }


class DescriptorLowerer:
    """
    Replacement seam for the verified EOS compiler.

    Required behavior:
        StageDescriptor6 -> EncodedQuintuplet
    """

    def lower_stage(self, descriptor: StageDescriptor6, stage_index: int, corner_name: str) -> EncodedQuintuplet:
        raise NotImplementedError(
            "DescriptorLowerer.lower_stage() must be filled from verified EOS decompilation"
        )

    def lower_bank(self, descriptors: Sequence[StageDescriptor6], corner_name: str) -> EncodedStageBank:
        if len(descriptors) != STAGE_COUNT:
            raise ValueError(f"expected {STAGE_COUNT} descriptors, got {len(descriptors)}")
        stages = tuple(
            self.lower_stage(desc, stage_index=i + 1, corner_name=corner_name)
            for i, desc in enumerate(descriptors)
        )
        return EncodedStageBank(stages)  # type: ignore[arg-type]


class FixtureLowerer(DescriptorLowerer):
    """
    Testing-only lowerer driven by explicit encoded fixture data.

    This lets the rest of the pipeline be tested before the real lowering code is
    wired in. It does not infer coefficients from descriptors.
    """

    def __init__(self, encoded_by_corner: Mapping[str, EncodedStageBank]) -> None:
        self.encoded_by_corner = dict(encoded_by_corner)

    def lower_stage(self, descriptor: StageDescriptor6, stage_index: int, corner_name: str) -> EncodedQuintuplet:
        try:
            return self.encoded_by_corner[corner_name].stages[stage_index - 1]
        except KeyError as exc:
            raise KeyError(f"fixture missing encoded corner {corner_name!r}") from exc


class EncodedDecoder:
    """Replacement seam for the verified EOS uint16 -> float c0..c4 decoder."""

    def decode_quintuplet(
        self,
        encoded: EncodedQuintuplet,
        host_sample_rate_hz: float = NATIVE_SAMPLE_RATE_HZ,
    ) -> Tuple[float, float, float, float, float]:
        raise NotImplementedError("EncodedDecoder.decode_quintuplet() must be filled from verified EOS logic")

    def decode_bank(
        self,
        bank: EncodedStageBank,
        host_sample_rate_hz: float = NATIVE_SAMPLE_RATE_HZ,
    ) -> List[Tuple[float, float, float, float, float]]:
        return [self.decode_quintuplet(stage, host_sample_rate_hz=host_sample_rate_hz) for stage in bank.stages]


class RossumRenderKernel:
    """Replacement seam for the verified proprietary serial biquad equation."""

    def process_sample(self, x: float, decoded_bank: Sequence[Tuple[float, float, float, float, float]]) -> float:
        raise NotImplementedError("RossumRenderKernel.process_sample() must be filled with the verified kernel equation")


class EosPipelineCompiler:
    """Strict four-stage EOS-style compiler orchestration."""

    def __init__(self, lowerer: DescriptorLowerer) -> None:
        self.lowerer = lowerer

    def compile_morph_designer_body(self, body: LogicalDescriptorBody) -> RuntimeEncodedSurface:
        m0_bank = self.lowerer.lower_bank(body.m0, corner_name="m0")
        m1_bank = self.lowerer.lower_bank(body.m1, corner_name="m1")
        return RuntimeEncodedSurface.from_morph_designer_pair(m0_bank, m1_bank)


def read_json_file(path: Path) -> Dict[str, Any]:
    resolved = resolve_path(path)
    with open(resolved, "r", encoding="utf-8") as f:
        data = json.load(f)
    if not isinstance(data, dict):
        raise ValueError(f"{resolved} must contain a JSON object")
    return data


def load_encoded_bank(matrix: Any, label: str) -> EncodedStageBank:
    if not isinstance(matrix, list):
        raise ValueError(f"{label}: expected list matrix")
    return EncodedStageBank.from_matrix(matrix)


def load_surface_from_fixture(path: Path) -> Tuple[str, RuntimeEncodedSurface, Dict[str, Any]]:
    data = read_json_file(path)
    name = str(data.get("name", "Untitled_EOS_Surface"))

    if "encoded_logical" in data:
        logical = data["encoded_logical"]
        if not isinstance(logical, Mapping):
            raise ValueError("encoded_logical must be an object with m0/m1")
        m0 = load_encoded_bank(logical.get("m0"), "encoded_logical.m0")
        m1 = load_encoded_bank(logical.get("m1"), "encoded_logical.m1")
        return name, RuntimeEncodedSurface.from_morph_designer_pair(m0, m1), data

    if "encoded_runtime" in data:
        runtime = data["encoded_runtime"]
        if not isinstance(runtime, Mapping):
            raise ValueError("encoded_runtime must be an object with four corners")
        corners = {name: load_encoded_bank(runtime.get(name), f"encoded_runtime.{name}") for name in RUNTIME_CORNER_ORDER}
        return name, RuntimeEncodedSurface.from_runtime_corners(corners), data

    raise ValueError("fixture must contain either encoded_logical or encoded_runtime")


def read_existing_revision(out_path: Path) -> int:
    resolved = resolve_path(out_path)
    try:
        if not resolved.exists():
            return 0
        with open(resolved, "r", encoding="utf-8") as f:
            data = json.load(f)
        return int(data.get("compiler", {}).get("revision", 0))
    except Exception:
        return 0


def make_hash_payload(payload: Mapping[str, Any]) -> Dict[str, Any]:
    copy = json.loads(stable_json_dumps(payload))
    compiler = copy.get("compiler")
    if isinstance(compiler, dict):
        compiler.pop("generated_at_utc", None)
        compiler.pop("revision", None)
        compiler.pop("content_hash", None)
    return copy


def make_cartridge_payload(
    *,
    name: str,
    surface: RuntimeEncodedSurface,
    revision: Optional[int],
    source: Mapping[str, Any],
    host_sample_rate_hz: float = NATIVE_SAMPLE_RATE_HZ,
) -> Dict[str, Any]:
    payload: Dict[str, Any] = {
        "schema": SCHEMA,
        "name": name,
        "engine_model": {
            "native_sample_rate_hz": NATIVE_SAMPLE_RATE_HZ,
            "host_sample_rate_hz": host_sample_rate_hz,
            "sqrt_compensation_level_inferred": sqrt_compensation_level(host_sample_rate_hz),
            "pipeline": [
                "descriptor_compiler",
                "encoded_uint16_bilinear_interpolator",
                "uint16_to_float_c_decoder",
                "proprietary_rossum_serial_6_stage_render_kernel",
            ],
            "interpolation_domain": "encoded_uint16_c_domain",
            "stage_count": STAGE_COUNT,
            "coeffs_per_stage": COEFFS_PER_STAGE,
            "q_authoring_note": (
                "Morph Designer authored two logical corners; Q rows are duplicated "
                "unless full runtime banks are supplied."
            ),
        },
        "compiler": {
            "name": COMPILER_NAME,
            "version": COMPILER_VERSION,
            "generated_at_utc": now_iso_utc(),
            "revision": revision,
        },
        "encoded_surface": {
            "corner_order": list(RUNTIME_CORNER_ORDER),
            "corners": surface.to_json_corners(),
        },
        "source": source,
    }
    payload["compiler"]["content_hash"] = short_hash(make_hash_payload(payload))
    return payload


def atomic_write_json(data: Mapping[str, Any], out_path: Path, keep_backup: bool = True) -> None:
    resolved = resolve_path(out_path)
    resolved.parent.mkdir(parents=True, exist_ok=True)

    tmp_path = resolved.with_name(resolved.name + ".tmp")
    backup_path = resolved.with_name(resolved.stem + ".last_good" + resolved.suffix)

    if keep_backup and resolved.exists():
        try:
            shutil.copy2(resolved, backup_path)
        except Exception:
            pass

    encoded = stable_json_dumps(data) + "\n"
    with open(tmp_path, "w", encoding="utf-8", newline="\n") as f:
        f.write(encoded)
        f.flush()
        os.fsync(f.fileno())

    os.replace(tmp_path, resolved)


def compile_fixture_to_file(
    fixture_path: Path,
    out_path: Path,
    *,
    host_sample_rate_hz: float = NATIVE_SAMPLE_RATE_HZ,
    keep_backup: bool = True,
) -> Dict[str, Any]:
    name, surface, fixture_data = load_surface_from_fixture(fixture_path)
    payload = make_cartridge_payload(
        name=name,
        surface=surface,
        revision=read_existing_revision(out_path) + 1,
        host_sample_rate_hz=host_sample_rate_hz,
        source={
            "type": "encoded_fixture",
            "path": str(resolve_path(fixture_path)),
            "fixture_keys": sorted(fixture_data.keys()),
        },
    )
    atomic_write_json(payload, out_path, keep_backup=keep_backup)
    return payload


EXAMPLE_FIXTURE = {
    "name": "example_encoded_fixture_replace_with_real_eos_output",
    "encoded_logical": {
        "m0": [[32768, 32768, 32768, 32768, 32768] for _ in range(STAGE_COUNT)],
        "m1": [[32768, 32768, 32768, 32768, 32768] for _ in range(STAGE_COUNT)],
    },
    "warning": "Example only. 32768 is not asserted to be EOS identity.",
}


def write_example_fixture(path: Path = DEFAULT_FIXTURE_PATH) -> Path:
    resolved = resolve_path(path)
    resolved.parent.mkdir(parents=True, exist_ok=True)
    with open(resolved, "w", encoding="utf-8", newline="\n") as f:
        f.write(stable_json_dumps(EXAMPLE_FIXTURE) + "\n")
    return resolved


def make_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="trench_eos_pipeline_compiler.py",
        description="EOS/TRENCH encoded-surface compiler harness.",
    )
    sub = parser.add_subparsers(dest="command")

    fixture_p = sub.add_parser("fixture", help="Compile/write from explicit encoded fixture JSON.")
    fixture_p.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE_PATH)
    fixture_p.add_argument("--out", type=Path, default=DEFAULT_OUT_PATH)
    fixture_p.add_argument("--host-sr", type=float, default=NATIVE_SAMPLE_RATE_HZ)
    fixture_p.add_argument("--no-backup", action="store_true")

    inspect_p = sub.add_parser("inspect", help="Inspect fixture output JSON without writing.")
    inspect_p.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE_PATH)
    inspect_p.add_argument("--host-sr", type=float, default=NATIVE_SAMPLE_RATE_HZ)

    interp_p = sub.add_parser("interp", help="Interpolate encoded fixture at Morph/Q and print encoded 6x5 bank.")
    interp_p.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE_PATH)
    interp_p.add_argument("--morph", type=float, default=0.5, help="0..1")
    interp_p.add_argument("--q", type=float, default=0.0, help="0..1")

    gui_p = sub.add_parser("gui", help="Open simple Tkinter fixture writer GUI.")
    gui_p.add_argument("--fixture", type=Path, default=DEFAULT_FIXTURE_PATH)
    gui_p.add_argument("--out", type=Path, default=DEFAULT_OUT_PATH)

    example_p = sub.add_parser("example-fixture", help="Write an example encoded fixture.")
    example_p.add_argument("--out", type=Path, default=DEFAULT_FIXTURE_PATH)

    return parser


def run_cli(argv: Optional[List[str]] = None) -> int:
    parser = make_arg_parser()
    args = parser.parse_args(argv)

    if args.command is None:
        return run_gui(DEFAULT_FIXTURE_PATH, DEFAULT_OUT_PATH)

    if args.command == "fixture":
        payload = compile_fixture_to_file(
            fixture_path=args.fixture,
            out_path=args.out,
            host_sample_rate_hz=args.host_sr,
            keep_backup=not args.no_backup,
        )
        print(f"Wrote {resolve_path(args.out)}")
        print(f"revision={payload['compiler']['revision']} hash={payload['compiler']['content_hash']}")
        return 0

    if args.command == "inspect":
        name, surface, _fixture_data = load_surface_from_fixture(args.fixture)
        payload = make_cartridge_payload(
            name=name,
            surface=surface,
            revision=None,
            host_sample_rate_hz=args.host_sr,
            source={"type": "encoded_fixture", "path": str(resolve_path(args.fixture))},
        )
        print(stable_json_dumps(payload))
        return 0

    if args.command == "interp":
        _name, surface, _fixture_data = load_surface_from_fixture(args.fixture)
        bank = surface.interpolate_encoded(morph=args.morph, q=args.q)
        print(stable_json_dumps({"morph": args.morph, "q": args.q, "interpolated_encoded_bank": bank.to_matrix()}))
        return 0

    if args.command == "gui":
        return run_gui(args.fixture, args.out)

    if args.command == "example-fixture":
        print(f"Wrote {write_example_fixture(args.out)}")
        return 0

    parser.error(f"Unknown command: {args.command}")
    return 2


def run_gui(fixture_path: Path, out_path: Path) -> int:
    try:
        import tkinter as tk
        from tkinter import filedialog, messagebox, ttk
    except Exception as exc:
        print(f"Tkinter unavailable: {exc}", file=sys.stderr)
        return 1

    class App:
        def __init__(self, root: "tk.Tk") -> None:
            self.root = root
            self.root.title("TRENCH EOS Pipeline Compiler")
            self.root.geometry("620x330")
            self.root.minsize(600, 300)
            self.root.configure(padx=18, pady=16)

            self.fixture_var = tk.StringVar(value=str(resolve_path(fixture_path)))
            self.out_var = tk.StringVar(value=str(resolve_path(out_path)))
            self.host_sr_var = tk.DoubleVar(value=NATIVE_SAMPLE_RATE_HZ)
            self.status_var = tk.StringVar(value="Ready. Load encoded fixture; no fake descriptor lowering is performed.")

            self._build(ttk, filedialog)

        def _build(self, ttk_module: Any, filedialog_module: Any) -> None:
            ttk_module.Label(self.root, text="TRENCH EOS Pipeline Compiler", font=("Segoe UI", 13, "bold")).pack(anchor="w")
            ttk_module.Label(
                self.root,
                text="Encoded uint16 c-domain surface writer. Lowerer/decoder/kernel must come from verified research.",
                foreground="gray",
                wraplength=580,
            ).pack(anchor="w", pady=(0, 14))

            self._path_row(ttk_module, filedialog_module, "Fixture", self.fixture_var, open_file=True)
            self._path_row(ttk_module, filedialog_module, "Output", self.out_var, open_file=False)

            sr_row = ttk_module.Frame(self.root)
            sr_row.pack(fill="x", pady=(6, 10))
            ttk_module.Label(sr_row, text="Host SR", width=10).pack(side="left")
            ttk_module.Entry(sr_row, textvariable=self.host_sr_var, width=14).pack(side="left")
            ttk_module.Label(sr_row, text=f"native = {NATIVE_SAMPLE_RATE_HZ} Hz", foreground="gray").pack(side="left", padx=(10, 0))

            btns = ttk_module.Frame(self.root)
            btns.pack(fill="x", pady=(8, 8))
            ttk_module.Button(btns, text="Validate Fixture", command=self.validate_fixture).pack(side="left")
            ttk_module.Button(btns, text="Compile & Save", command=self.compile_save).pack(side="left", padx=(8, 0), fill="x", expand=True)
            ttk_module.Button(btns, text="Print Interp 50%", command=self.print_interp).pack(side="left", padx=(8, 0))

            ttk_module.Separator(self.root).pack(fill="x", pady=(10, 8))
            ttk_module.Label(self.root, textvariable=self.status_var, foreground="gray", wraplength=580, justify="left").pack(anchor="w")

        def _path_row(self, ttk_module: Any, filedialog_module: Any, label: str, var: Any, *, open_file: bool) -> None:
            row = ttk_module.Frame(self.root)
            row.pack(fill="x", pady=(0, 8))
            ttk_module.Label(row, text=label, width=10).pack(side="left")
            ttk_module.Entry(row, textvariable=var).pack(side="left", fill="x", expand=True)

            def browse() -> None:
                if open_file:
                    selected = filedialog_module.askopenfilename(filetypes=[("JSON", "*.json"), ("All files", "*.*")])
                else:
                    selected = filedialog_module.asksaveasfilename(defaultextension=".json", filetypes=[("JSON", "*.json"), ("All files", "*.*")])
                if selected:
                    var.set(selected)

            ttk_module.Button(row, text="Browse", command=browse).pack(side="left", padx=(8, 0))

        def _fixture_path(self) -> Path:
            return Path(self.fixture_var.get().strip())

        def _out_path(self) -> Path:
            return Path(self.out_var.get().strip())

        def validate_fixture(self) -> None:
            try:
                name, surface, _ = load_surface_from_fixture(self._fixture_path())
                self.status_var.set(
                    f"Valid fixture: {name}\n"
                    f"Corners: {', '.join(surface.to_json_corners().keys())}\n"
                    f"Q rows duplicated: {surface.m0_q0 == surface.m0_q1 and surface.m1_q0 == surface.m1_q1}"
                )
            except Exception as exc:
                messagebox.showerror("Fixture Error", str(exc))
                self.status_var.set(f"Fixture invalid: {exc}")

        def compile_save(self) -> None:
            try:
                payload = compile_fixture_to_file(
                    fixture_path=self._fixture_path(),
                    out_path=self._out_path(),
                    host_sample_rate_hz=float(self.host_sr_var.get()),
                    keep_backup=True,
                )
                self.status_var.set(
                    f"Wrote revision {payload['compiler']['revision']} -> {resolve_path(self._out_path())}\n"
                    f"hash={payload['compiler']['content_hash']}"
                )
            except Exception as exc:
                messagebox.showerror("Compile Error", str(exc))
                self.status_var.set(f"Compile failed: {exc}")

        def print_interp(self) -> None:
            try:
                _name, surface, _ = load_surface_from_fixture(self._fixture_path())
                bank = surface.interpolate_encoded(morph=0.5, q=0.0)
                print(stable_json_dumps({"morph": 0.5, "q": 0.0, "bank": bank.to_matrix()}))
                self.status_var.set("Printed encoded 50% Morph / Q0 interpolation bank to console.")
            except Exception as exc:
                messagebox.showerror("Interpolation Error", str(exc))
                self.status_var.set(f"Interpolation failed: {exc}")

    root = tk.Tk()
    style = ttk.Style()
    if "clam" in style.theme_names():
        style.theme_use("clam")
    App(root)
    root.mainloop()
    return 0


if __name__ == "__main__":
    raise SystemExit(run_cli())
