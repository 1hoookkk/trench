"""Shipping finalists — tables-driven candidate generation for the 4 flagship bodies.

Generates exactly 4 editable candidates per flagship body using:
- SHIPPING.md body law
- docs/sonic_tables/tables.json priors
- existing direct StageParams pipeline

This is not the old bulk Pareto dump. It emits a clean finalist set only.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, replace
from pathlib import Path

from pyruntime.body import Body
from pyruntime.forge_generator import cascade_peak_db
from pyruntime.forge_shipping import (
    ROLE_ANCHOR,
    ROLE_BOUNDARY,
    ROLE_CHARACTER,
    ShippingBodyDef,
    StageTarget,
    build_body_from_def,
    eval_shipping_body,
)

ROOT = Path(__file__).parent.parent
TABLES_PATH = ROOT / "docs" / "sonic_tables" / "tables.json"
VAULT_DIR = ROOT / "vault"
FINAL_DIR = VAULT_DIR / "_shipping_finalists"


@dataclass(frozen=True)
class CandidatePackage:
    body_key: str
    candidate_key: str
    body_def: ShippingBodyDef
    table_refs: list[str]
    notes: str


def load_tables() -> dict:
    return json.loads(TABLES_PATH.read_text())


def vowel(tables: dict, key: str) -> dict:
    return tables["vowels"][key]


def nasal(tables: dict, key: str) -> dict:
    return tables["nasals"][key]


def landmark(tables: dict, name: str) -> dict:
    for entry in tables["landmarks"]["entries"]:
        if entry["name"] == name:
            return entry
    raise KeyError(name)


def consonant(tables: dict, key: str) -> dict:
    return tables["consonants"][key]


def instrument(tables: dict, key: str) -> dict:
    return tables["instruments"][key]


def moog(tables: dict, key: str) -> dict:
    return tables["moog"][key]


def singer(tables: dict) -> dict:
    return tables["singer"]["singers_formant"]


def bell(tables: dict, index: int) -> dict:
    return tables["measured_bells"][index]


def candidate_definitions() -> list[CandidatePackage]:
    tables = load_tables()

    oo = vowel(tables, "oo")
    ah = vowel(tables, "ah")
    ee = vowel(tables, "ee")
    schwa = vowel(tables, "schwa")
    nasal_m = nasal(tables, "nasal_m")
    chest = landmark(tables, "chest_resonance")
    room_medium = landmark(tables, "room_mode_medium")
    room_small = landmark(tables, "room_mode_small")
    telephone_low = landmark(tables, "telephone_band_low")
    telephone_high = landmark(tables, "telephone_band_high")
    presence = landmark(tables, "presence_peak")
    sibilance = landmark(tables, "sibilance_band")
    air = landmark(tables, "air_band")
    s = consonant(tables, "s")
    sh = consonant(tables, "sh")
    trumpet = instrument(tables, "trumpet")
    oboe = instrument(tables, "oboe")
    ladder = moog(tables, "moog_uncompensated")
    singer_ring = singer(tables)
    bell_a = bell(tables, 0)
    bell_b = bell(tables, 3)

    return [
        CandidatePackage(
            "speaker_knockerz",
            "cand_01_vault_choke",
            ShippingBodyDef(
                name="Speaker Knockerz C1",
                key="speaker_knockerz",
                m0_q0=[
                    StageTarget(45.0, 0.979, ROLE_ANCHOR),
                    StageTarget(float(room_medium["freq_hz"]), 0.968, ROLE_CHARACTER),
                    StageTarget(float(chest["freq_hz"]), 0.970, ROLE_CHARACTER, "interior"),
                    StageTarget(250.0, 0.966, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(400.0, 0.964, ROLE_CHARACTER, "interior"),
                    StageTarget(float(air["freq_hz"]), 0.942, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(45.0, 0.979, ROLE_ANCHOR),
                    StageTarget(float(chest["freq_hz"]), 0.979, ROLE_CHARACTER),
                    StageTarget(220.0, 0.976, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(400.0, 0.983, ROLE_CHARACTER, "interior"),
                    StageTarget(1200.0, 0.982, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(float(sibilance["freq_hz"]), 0.946, ROLE_BOUNDARY),
                ],
                q_compression=0.52,
                q_center_hz=420.0,
                q_radius_mod=0.004,
                min_radius=0.958,
            ),
            [
                "landmarks.room_mode_medium",
                "landmarks.chest_resonance",
                "landmarks.air_band",
            ],
            "Sub stays pinned while the low-mid throat tightens into choke and cone cry.",
        ),
        CandidatePackage(
            "speaker_knockerz",
            "cand_02_cardboard_rip",
            ShippingBodyDef(
                name="Speaker Knockerz C2",
                key="speaker_knockerz",
                m0_q0=[
                    StageTarget(45.0, 0.978, ROLE_ANCHOR),
                    StageTarget(float(room_small["freq_hz"]), 0.969, ROLE_CHARACTER),
                    StageTarget(float(telephone_low["freq_hz"]), 0.966, ROLE_CHARACTER, "interior"),
                    StageTarget(500.0, 0.963, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(float(ladder["fc_default"]), 0.960, ROLE_CHARACTER, "interior"),
                    StageTarget(float(air["freq_hz"]), 0.940, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(45.0, 0.978, ROLE_ANCHOR),
                    StageTarget(180.0, 0.977, ROLE_CHARACTER),
                    StageTarget(320.0, 0.980, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(780.0, 0.983, ROLE_CHARACTER, "interior"),
                    StageTarget(2500.0, 0.981, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(float(sibilance["freq_hz"]), 0.944, ROLE_BOUNDARY),
                ],
                q_compression=0.48,
                q_center_hz=380.0,
                q_radius_mod=0.005,
                min_radius=0.956,
            ),
            [
                "landmarks.room_mode_small",
                "landmarks.telephone_band_low",
                "moog.moog_uncompensated",
            ],
            "More paper-cone aggression and more obvious phone-speaker translation damage.",
        ),
        CandidatePackage(
            "speaker_knockerz",
            "cand_03_rattle_howl",
            ShippingBodyDef(
                name="Speaker Knockerz C3",
                key="speaker_knockerz",
                m0_q0=[
                    StageTarget(45.0, 0.980, ROLE_ANCHOR),
                    StageTarget(95.0, 0.968, ROLE_CHARACTER),
                    StageTarget(float(chest["freq_hz"]), 0.969, ROLE_CHARACTER, "interior"),
                    StageTarget(300.0, 0.965, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(820.0, 0.962, ROLE_CHARACTER, "interior"),
                    StageTarget(float(air["freq_hz"]), 0.939, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(45.0, 0.980, ROLE_ANCHOR),
                    StageTarget(150.0, 0.978, ROLE_CHARACTER),
                    StageTarget(420.0, 0.981, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(950.0, 0.984, ROLE_CHARACTER, "interior"),
                    StageTarget(1800.0, 0.982, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(9000.0, 0.944, ROLE_BOUNDARY),
                ],
                q_compression=0.56,
                q_center_hz=460.0,
                q_radius_mod=0.004,
                min_radius=0.958,
            ),
            [
                "landmarks.chest_resonance",
                "landmarks.air_band",
            ],
            "A little less tear, a little more mechanical rattle and cone howl.",
        ),
        CandidatePackage(
            "speaker_knockerz",
            "cand_04_fracture_hiss",
            ShippingBodyDef(
                name="Speaker Knockerz C4",
                key="speaker_knockerz",
                m0_q0=[
                    StageTarget(45.0, 0.978, ROLE_ANCHOR),
                    StageTarget(float(room_medium["freq_hz"]), 0.967, ROLE_CHARACTER),
                    StageTarget(170.0, 0.968, ROLE_CHARACTER, "interior"),
                    StageTarget(400.0, 0.964, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(1200.0, 0.960, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(float(sibilance["freq_hz"]), 0.942, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(45.0, 0.978, ROLE_ANCHOR),
                    StageTarget(125.0, 0.976, ROLE_CHARACTER),
                    StageTarget(260.0, 0.979, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(640.0, 0.983, ROLE_CHARACTER, "interior"),
                    StageTarget(2400.0, 0.982, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(float(air["freq_hz"]), 0.946, ROLE_BOUNDARY),
                ],
                q_compression=0.54,
                q_center_hz=360.0,
                q_radius_mod=0.006,
                min_radius=0.956,
            ),
            [
                "landmarks.room_mode_medium",
                "landmarks.sibilance_band",
                "landmarks.air_band",
            ],
            "Most hostile upper fracture while still keeping the sub pinned.",
        ),
        CandidatePackage(
            "aluminum_siding",
            "cand_01_sheen_fold",
            ShippingBodyDef(
                name="Aluminum Siding C1",
                key="aluminum_siding",
                m0_q0=[
                    StageTarget(300.0, 0.960, ROLE_BOUNDARY),
                    StageTarget(float(presence["freq_hz"]), 0.972, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(float(sibilance["freq_hz"]), 0.975, ROLE_CHARACTER),
                    StageTarget(float(air["freq_hz"]), 0.971, ROLE_CHARACTER),
                    StageTarget(float(sh["formants"][3]["freq"]), 0.969, ROLE_CHARACTER, "interior"),
                    StageTarget(15000.0, 0.958, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(300.0, 0.960, ROLE_BOUNDARY),
                    StageTarget(5000.0, 0.974, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(7200.0, 0.977, ROLE_CHARACTER),
                    StageTarget(11000.0, 0.974, ROLE_CHARACTER),
                    StageTarget(15000.0, 0.970, ROLE_CHARACTER, "interior"),
                    StageTarget(17500.0, 0.954, ROLE_BOUNDARY),
                ],
                q_compression=0.38,
                q_center_hz=8500.0,
                q_radius_mod=0.004,
                min_radius=0.960,
            ),
            [
                "landmarks.presence_peak",
                "landmarks.sibilance_band",
                "landmarks.air_band",
                "consonants.sh",
            ],
            "Starts glossy, then folds the sibilance band inward before it tears upward.",
        ),
        CandidatePackage(
            "aluminum_siding",
            "cand_02_glass_whistle",
            ShippingBodyDef(
                name="Aluminum Siding C2",
                key="aluminum_siding",
                m0_q0=[
                    StageTarget(263.6, 0.898, ROLE_BOUNDARY, "unit_circle"),
                    StageTarget(1116.3, 0.900, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(7593.4, 0.959, ROLE_CHARACTER, "interior"),
                    StageTarget(8815.5, 0.915, ROLE_CHARACTER, "interior"),
                    StageTarget(14420.2, 0.904, ROLE_CHARACTER, "interior"),
                    StageTarget(18936.8, 0.969, ROLE_BOUNDARY, "interior"),
                ],
                m100_q0=[
                    StageTarget(280.2, 0.947, ROLE_BOUNDARY, "unit_circle"),
                    StageTarget(970.6, 0.927, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(7655.7, 0.975, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(7002.3, 0.960, ROLE_CHARACTER, "interior"),
                    StageTarget(12378.5, 0.882, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(18016.3, 0.929, ROLE_BOUNDARY, "unit_circle"),
                ],
                q_compression=0.32015786686717856,
                q_center_hz=9453.638395144453,
                q_radius_mod=0.00012293458881456164,
                min_radius=0.961,
            ),
            [
                "consonants.s",
                "landmarks.sibilance_band",
                "landmarks.air_band",
            ],
            "More glass and whistle, less silk.",
        ),
        CandidatePackage(
            "aluminum_siding",
            "cand_03_sibilant_swallow",
            ShippingBodyDef(
                name="Aluminum Siding C3",
                key="aluminum_siding",
                m0_q0=[
                    StageTarget(300.0, 0.960, ROLE_BOUNDARY),
                    StageTarget(5000.0, 0.973, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(8000.0, 0.974, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(9600.0, 0.970, ROLE_CHARACTER),
                    StageTarget(float(air["freq_hz"]), 0.968, ROLE_CHARACTER),
                    StageTarget(16500.0, 0.956, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(300.0, 0.960, ROLE_BOUNDARY),
                    StageTarget(5400.0, 0.975, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(8200.0, 0.977, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(12500.0, 0.974, ROLE_CHARACTER),
                    StageTarget(17000.0, 0.971, ROLE_CHARACTER),
                    StageTarget(18500.0, 0.952, ROLE_BOUNDARY),
                ],
                q_compression=0.41,
                q_center_hz=8800.0,
                q_radius_mod=0.004,
                min_radius=0.960,
            ),
            [
                "landmarks.sibilance_band",
                "landmarks.air_band",
            ],
            "More obvious S/T self-swallowing fold, less pop-air polish.",
        ),
        CandidatePackage(
            "aluminum_siding",
            "cand_04_foil_shatter",
            ShippingBodyDef(
                name="Aluminum Siding C4",
                key="aluminum_siding",
                m0_q0=[
                    StageTarget(300.0, 0.960, ROLE_BOUNDARY),
                    StageTarget(float(presence["freq_hz"]), 0.971, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(7000.0, 0.974, ROLE_CHARACTER),
                    StageTarget(10000.0, 0.972, ROLE_CHARACTER),
                    StageTarget(13000.0, 0.970, ROLE_CHARACTER, "interior"),
                    StageTarget(17000.0, 0.955, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(300.0, 0.960, ROLE_BOUNDARY),
                    StageTarget(4700.0, 0.974, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(9000.0, 0.977, ROLE_CHARACTER),
                    StageTarget(14000.0, 0.974, ROLE_CHARACTER),
                    StageTarget(18200.0, 0.972, ROLE_CHARACTER, "interior"),
                    StageTarget(19400.0, 0.952, ROLE_BOUNDARY),
                ],
                q_compression=0.46,
                q_center_hz=9600.0,
                q_radius_mod=0.005,
                min_radius=0.962,
            ),
            [
                "landmarks.presence_peak",
                "landmarks.air_band",
            ],
            "Most metallic foil-tear and shatter-point candidate.",
        ),
        CandidatePackage(
            "small_talk",
            "cand_01_open_bite",
            ShippingBodyDef(
                name="Small Talk Ah-Ee C1",
                key="small_talk",
                m0_q0=[
                    StageTarget(float(ah["f1"]), 0.972, ROLE_CHARACTER),
                    StageTarget(float(ah["f2"]), 0.971, ROLE_CHARACTER),
                    StageTarget(float(schwa["f3"]), 0.955, ROLE_BOUNDARY, "interior"),
                    StageTarget(180.0, 0.950, ROLE_BOUNDARY, "bandpass"),
                    StageTarget(float(telephone_high["freq_hz"]), 0.945, ROLE_BOUNDARY),
                    StageTarget(4800.0, 0.938, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(float(ee["f1"]), 0.973, ROLE_CHARACTER),
                    StageTarget(float(ee["f2"]), 0.972, ROLE_CHARACTER),
                    StageTarget(float(singer_ring["f4"]), 0.958, ROLE_BOUNDARY, "interior"),
                    StageTarget(float(nasal_m["f2"]), 0.950, ROLE_BOUNDARY, "bandpass"),
                    StageTarget(float(telephone_high["freq_hz"]), 0.944, ROLE_BOUNDARY),
                    StageTarget(5200.0, 0.936, ROLE_BOUNDARY),
                ],
                q_compression=0.20,
                q_center_hz=1350.0,
                q_radius_mod=0.0,
                min_radius=0.940,
            ),
            [
                "vowels.ah",
                "vowels.ee",
                "vowels.schwa",
                "nasals.nasal_m",
                "singer.singers_formant",
            ],
            "Direct Ah to Ee with a nasal undercarriage and upper bite ring.",
        ),
        CandidatePackage(
            "small_talk",
            "cand_02_yawn_to_bite",
            ShippingBodyDef(
                name="Small Talk Ah-Ee C2",
                key="small_talk",
                m0_q0=[
                    StageTarget(float(oo["f1"]), 0.971, ROLE_CHARACTER),
                    StageTarget(float(oo["f2"]), 0.970, ROLE_CHARACTER),
                    StageTarget(float(oo["f3"]), 0.954, ROLE_BOUNDARY, "interior"),
                    StageTarget(float(nasal_m["f1"]), 0.948, ROLE_BOUNDARY, "bandpass"),
                    StageTarget(float(telephone_high["freq_hz"]), 0.944, ROLE_BOUNDARY),
                    StageTarget(4700.0, 0.936, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(float(ee["f1"]), 0.973, ROLE_CHARACTER),
                    StageTarget(float(ee["f2"]), 0.972, ROLE_CHARACTER),
                    StageTarget(float(ee["f3"]), 0.958, ROLE_BOUNDARY, "interior"),
                    StageTarget(float(singer_ring["f3"]), 0.950, ROLE_BOUNDARY),
                    StageTarget(float(telephone_high["freq_hz"]), 0.944, ROLE_BOUNDARY),
                    StageTarget(5200.0, 0.936, ROLE_BOUNDARY),
                ],
                q_compression=0.20,
                q_center_hz=1180.0,
                q_radius_mod=0.0,
                min_radius=0.940,
            ),
            [
                "vowels.oo",
                "vowels.ee",
                "nasals.nasal_m",
                "singer.singers_formant",
            ],
            "Bigger cavity at the floor, tighter and more synthetic bite at the roof.",
        ),
        CandidatePackage(
            "small_talk",
            "cand_03_hum_gag",
            ShippingBodyDef(
                name="Small Talk Ah-Ee C3",
                key="small_talk",
                m0_q0=[
                    StageTarget(float(nasal_m["f1"]), 0.966, ROLE_CHARACTER),
                    StageTarget(float(nasal_m["f2"]), 0.965, ROLE_CHARACTER),
                    StageTarget(float(ah["f1"]), 0.952, ROLE_BOUNDARY),
                    StageTarget(1700.0, 0.948, ROLE_BOUNDARY, "interior"),
                    StageTarget(300.0, 0.944, ROLE_BOUNDARY, "bandpass"),
                    StageTarget(float(telephone_high["freq_hz"]), 0.936, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(420.0, 0.968, ROLE_CHARACTER),
                    StageTarget(float(ee["f2"]), 0.968, ROLE_CHARACTER),
                    StageTarget(float(singer_ring["f3"]), 0.952, ROLE_BOUNDARY),
                    StageTarget(float(singer_ring["f5"]), 0.948, ROLE_BOUNDARY, "interior"),
                    StageTarget(float(nasal_m["f2"]), 0.944, ROLE_BOUNDARY, "bandpass"),
                    StageTarget(float(telephone_high["freq_hz"]), 0.936, ROLE_BOUNDARY),
                ],
                q_compression=0.12,
                q_center_hz=1250.0,
                q_radius_mod=0.0,
                min_radius=0.940,
            ),
            [
                "nasals.nasal_m",
                "vowels.ah",
                "vowels.ee",
                "singer.singers_formant",
            ],
            "Most obvious hum-to-gag candidate, with stronger low-formant choke.",
        ),
        CandidatePackage(
            "small_talk",
            "cand_04_trumpet_tube",
            ShippingBodyDef(
                name="Small Talk Ah-Ee C4",
                key="small_talk",
                m0_q0=[
                    StageTarget(float(ah["f1"]), 0.968, ROLE_CHARACTER),
                    StageTarget(float(ah["f2"]), 0.967, ROLE_CHARACTER),
                    StageTarget(float(trumpet["peaks"][0]["freq"]), 0.952, ROLE_BOUNDARY),
                    StageTarget(float(oboe["peaks"][1]["freq"]), 0.949, ROLE_BOUNDARY),
                    StageTarget(180.0, 0.944, ROLE_BOUNDARY, "bandpass"),
                    StageTarget(float(telephone_high["freq_hz"]), 0.936, ROLE_BOUNDARY),
                ],
                m100_q0=[
                    StageTarget(float(ee["f1"]), 0.969, ROLE_CHARACTER),
                    StageTarget(float(ee["f2"]), 0.969, ROLE_CHARACTER),
                    StageTarget(float(trumpet["peaks"][1]["freq"]), 0.953, ROLE_BOUNDARY),
                    StageTarget(float(oboe["peaks"][0]["freq"]), 0.949, ROLE_BOUNDARY),
                    StageTarget(float(singer_ring["f3"]), 0.944, ROLE_BOUNDARY, "bandpass"),
                    StageTarget(float(telephone_high["freq_hz"]), 0.936, ROLE_BOUNDARY),
                ],
                q_compression=0.10,
                q_center_hz=1420.0,
                q_radius_mod=0.0,
                min_radius=0.940,
            ),
            [
                "vowels.ah",
                "vowels.ee",
                "instruments.trumpet",
                "instruments.oboe",
            ],
            "More synthetic tube pressure and brassy bite than the others.",
        ),
        CandidatePackage(
            "cul_de_sac",
            "cand_01_pipe_to_shards",
            ShippingBodyDef(
                name="Cul-De-Sac C1",
                key="cul_de_sac",
                m0_q0=[
                    StageTarget(100.0, 0.975, ROLE_ANCHOR),
                    StageTarget(260.0, 0.969, ROLE_CHARACTER),
                    StageTarget(390.0, 0.967, ROLE_CHARACTER),
                    StageTarget(520.0, 0.966, ROLE_CHARACTER),
                    StageTarget(700.0, 0.964, ROLE_CHARACTER, "interior"),
                    StageTarget(900.0, 0.962, ROLE_CHARACTER),
                ],
                m100_q0=[
                    StageTarget(100.0, 0.975, ROLE_ANCHOR),
                    StageTarget(float(bell_a["partials"][2]), 0.983, ROLE_CHARACTER),
                    StageTarget(float(telephone_high["freq_hz"]), 0.982, ROLE_CHARACTER),
                    StageTarget(float(bell_b["partials"][4]), 0.980, ROLE_CHARACTER),
                    StageTarget(float(sibilance["freq_hz"]), 0.978, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(float(air["freq_hz"]), 0.976, ROLE_CHARACTER),
                ],
                q_compression=0.60,
                q_center_hz=2100.0,
                q_radius_mod=0.004,
                min_radius=0.958,
            ),
            [
                "measured_bells.0",
                "measured_bells.3",
                "landmarks.telephone_band_high",
                "landmarks.sibilance_band",
            ],
            "Pipe cluster fractures into metallic shards with an audible null seam.",
        ),
        CandidatePackage(
            "cul_de_sac",
            "cand_02_null_bridge",
            ShippingBodyDef(
                name="Cul-De-Sac C2",
                key="cul_de_sac",
                m0_q0=[
                    StageTarget(100.0, 0.975, ROLE_ANCHOR),
                    StageTarget(300.0, 0.968, ROLE_CHARACTER),
                    StageTarget(420.0, 0.966, ROLE_CHARACTER),
                    StageTarget(540.0, 0.965, ROLE_CHARACTER),
                    StageTarget(700.0, 0.963, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(820.0, 0.961, ROLE_CHARACTER),
                ],
                m100_q0=[
                    StageTarget(100.0, 0.975, ROLE_ANCHOR),
                    StageTarget(1500.0, 0.984, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(3400.0, 0.983, ROLE_CHARACTER),
                    StageTarget(6200.0, 0.981, ROLE_CHARACTER),
                    StageTarget(10200.0, 0.980, ROLE_CHARACTER),
                    StageTarget(14500.0, 0.978, ROLE_CHARACTER),
                ],
                q_compression=0.63,
                q_center_hz=1800.0,
                q_radius_mod=0.005,
                min_radius=0.958,
            ),
            [
                "landmarks.telephone_band_high",
                "landmarks.sibilance_band",
                "landmarks.air_band",
            ],
            "Most obvious mid-morph null and the cleanest fracture jump.",
        ),
        CandidatePackage(
            "cul_de_sac",
            "cand_03_comb_dust",
            ShippingBodyDef(
                name="Cul-De-Sac C3",
                key="cul_de_sac",
                m0_q0=[
                    StageTarget(73.1, 0.848, ROLE_ANCHOR, "bandpass"),
                    StageTarget(172.8, 0.802, ROLE_ANCHOR, "bandpass"),
                    StageTarget(325.0, 0.885, ROLE_CHARACTER, "allpole"),
                    StageTarget(455.0, 0.873, ROLE_CHARACTER, "interior"),
                    StageTarget(697.5, 0.882, ROLE_CHARACTER, "unit_circle"),
                    StageTarget(3358.3, 0.940, ROLE_CHARACTER, "interior"),
                ],
                m100_q0=[
                    StageTarget(73.1, 0.851, ROLE_ANCHOR, "bandpass"),
                    StageTarget(172.8, 0.860, ROLE_ANCHOR, "allpole"),
                    StageTarget(1829.2, 0.939, ROLE_CHARACTER, "allpole"),
                    StageTarget(3297.8, 0.927, ROLE_CHARACTER, "interior"),
                    StageTarget(10054.2, 0.915, ROLE_CHARACTER, "allpole"),
                    StageTarget(14372.6, 0.950, ROLE_CHARACTER, "unit_circle"),
                ],
                q_compression=0.009512738463790183,
                q_center_hz=606.4594577101592,
                q_radius_mod=0.001099486378714862,
                min_radius=0.958,
            ),
            [
                "measured_bells.0",
                "measured_bells.3",
            ],
            "Less pipe, more comb dust and crystalline high-end separation.",
        ),
        CandidatePackage(
            "cul_de_sac",
            "cand_04_bulge_collapse",
            ShippingBodyDef(
                name="Cul-De-Sac C4",
                key="cul_de_sac",
                m0_q0=[
                    StageTarget(100.0, 0.975, ROLE_ANCHOR),
                    StageTarget(320.0, 0.969, ROLE_CHARACTER),
                    StageTarget(430.0, 0.967, ROLE_CHARACTER),
                    StageTarget(520.0, 0.966, ROLE_CHARACTER),
                    StageTarget(680.0, 0.964, ROLE_CHARACTER),
                    StageTarget(3000.0, 0.962, ROLE_CHARACTER, "unit_circle"),
                ],
                m100_q0=[
                    StageTarget(100.0, 0.975, ROLE_ANCHOR),
                    StageTarget(1700.0, 0.984, ROLE_CHARACTER),
                    StageTarget(3600.0, 0.983, ROLE_CHARACTER),
                    StageTarget(6000.0, 0.981, ROLE_CHARACTER),
                    StageTarget(11000.0, 0.980, ROLE_CHARACTER),
                    StageTarget(15000.0, 0.978, ROLE_CHARACTER, "unit_circle"),
                ],
                q_compression=0.62,
                q_center_hz=2000.0,
                q_radius_mod=0.005,
                min_radius=0.958,
            ),
            [
                "landmarks.telephone_band_high",
                "landmarks.air_band",
            ],
            "Most pronounced 3k bulge before the structure caves into shards.",
        ),
    ]


def _normalize_stage(body_key: str, stage: StageTarget) -> StageTarget:
    if body_key == "speaker_knockerz":
        zero_type = stage.zero_type
        if stage.freq_hz <= 300.0 and stage.role != ROLE_ANCHOR:
            zero_type = "bandpass"
        if stage.role == ROLE_ANCHOR:
            radius = min(stage.radius, 0.955)
        elif stage.freq_hz < 1500.0:
            radius = min(stage.radius, 0.955)
        else:
            radius = min(stage.radius, 0.945)
        return replace(stage, radius=radius, zero_type=zero_type)

    if body_key == "cul_de_sac":
        zero_type = stage.zero_type
        if stage.freq_hz <= 900.0 and stage.role != ROLE_ANCHOR:
            zero_type = "bandpass"
        if stage.role == ROLE_ANCHOR:
            radius = min(stage.radius, 0.952)
        elif stage.freq_hz < 1200.0:
            radius = min(stage.radius, 0.950)
        else:
            radius = min(stage.radius, 0.965)
        return replace(stage, radius=radius, zero_type=zero_type)

    return stage


def _compile_body_def(pkg: CandidatePackage) -> ShippingBodyDef:
    if pkg.body_key in {"speaker_knockerz"}:
        return replace(
            pkg.body_def,
            m0_q0=[_normalize_stage(pkg.body_key, stage) for stage in pkg.body_def.m0_q0],
            m100_q0=[_normalize_stage(pkg.body_key, stage) for stage in pkg.body_def.m100_q0],
            q_radius_mod=min(pkg.body_def.q_radius_mod, 0.002),
            target_db=24.0,
            min_radius=min(pkg.body_def.min_radius, 0.950),
        )
    return pkg.body_def


def compile_candidate(pkg: CandidatePackage) -> tuple[Body, dict]:
    compiled_def = _compile_body_def(pkg)
    body = build_body_from_def(compiled_def)
    body = Body(name=compiled_def.name, corners=body.corners, boost=body.boost)
    metrics = eval_shipping_body(body, pkg.body_key)
    metrics["peaks"] = cascade_peak_db(body)
    metrics["table_refs"] = pkg.table_refs
    metrics["notes"] = pkg.notes
    return body, metrics


def write_finalists(packages: list[CandidatePackage]) -> Path:
    FINAL_DIR.mkdir(parents=True, exist_ok=True)
    manifest: dict[str, list[dict]] = {}

    grouped: dict[str, list[CandidatePackage]] = {}
    for pkg in packages:
        grouped.setdefault(pkg.body_key, []).append(pkg)

    for body_key, body_packages in grouped.items():
        body_dir = FINAL_DIR / body_key
        body_dir.mkdir(parents=True, exist_ok=True)
        manifest[body_key] = []

        for pkg in body_packages:
            body, metrics = compile_candidate(pkg)
            out_name = f"{pkg.body_key}__{pkg.candidate_key}.json"
            out_path = body_dir / out_name
            out_path.write_text(body.to_compiled_json(provenance="shipping-finalists"))
            manifest[body_key].append(
                {
                    "candidate": pkg.candidate_key,
                    "name": pkg.body_def.name,
                    "path": str(out_path),
                    "metrics": metrics,
                    "table_refs": pkg.table_refs,
                    "notes": pkg.notes,
                }
            )

    (FINAL_DIR / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return FINAL_DIR


def main() -> None:
    out_dir = write_finalists(candidate_definitions())
    print(out_dir)


if __name__ == "__main__":
    main()
