use trench_core::cartridge::Cartridge;

#[test]
fn existing_cartridge_without_new_blocks_loads_with_safe_defaults() {
    let json = std::fs::read_to_string("../trench-juce/plugin/assets/cartridges/aluminum_siding.json")
        .expect("read fixture");
    let cart = Cartridge::from_json(&json).expect("parse");
    assert_eq!(cart.drive.input_gain_db, 0.0);
    assert_eq!(cart.drive.model, "mackie_1202");
    assert!(cart.spatial_profile.is_none());
    assert!(cart.mod_fn.is_none());
}

#[test]
fn cartridge_with_drive_mod_fn_spatial_parses() {
    // Minimum keyframes stub that passes the existing compiled-v1 parser.
    let json = r#"{
      "format": "compiled-v1",
      "name": "test",
      "sampleRate": 48000,
      "keyframes": [
        {"label": "M0_Q0", "stages": [
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}
        ]},
        {"label": "M100_Q0", "stages": [
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}
        ]},
        {"label": "M0_Q100", "stages": [
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}
        ]},
        {"label": "M100_Q100", "stages": [
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},
          {"c0":1,"c1":0,"c2":0,"c3":0,"c4":0},{"c0":1,"c1":0,"c2":0,"c3":0,"c4":0}
        ]}
      ],
      "drive": { "input_gain_dB": 9.0, "model": "mackie_1202" },
      "mod_fn": {
        "segments": [{ "level": 1.0, "time_ms": 50, "shape": "exp" }],
        "key-sync": 1, "tempo-sync": 1
      }
    }"#;
    let cart = Cartridge::from_json(json).expect("parse");
    assert_eq!(cart.drive.input_gain_db, 9.0);
    let mf = cart.mod_fn.as_ref().expect("mod_fn present");
    assert_eq!(mf.segments[0].level, 1.0);
    assert!(mf.key_sync());
    assert!(mf.tempo_sync());
}
