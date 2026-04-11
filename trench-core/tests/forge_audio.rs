/// Integration test: load a forge body, process impulse and noise, confirm non-silent output.
use trench_core::{Cartridge, Cascade};

fn load_acid_squelch() -> Cartridge {
    let json = std::fs::read_to_string("../cartridges/Acid_Squelch_77.json")
        .expect("../cartridges/Acid_Squelch_77.json must exist — run forge_generator first");
    Cartridge::from_json(&json).expect("failed to parse forge body")
}

#[test]
fn forge_body_impulse_produces_output() {
    let cart = load_acid_squelch();
    let mut cascade = Cascade::new();

    // Set coefficients at morph=0.5, q=0.0
    let coeffs = cart.interpolate(0.5, 0.0);

    // Feed an impulse: 1.0 followed by silence
    let mut buf = [0.0f32; 256];
    buf[0] = 1.0;
    for chunk in buf.chunks_mut(trench_core::BLOCK_SIZE) {
        cascade.set_targets(&coeffs, chunk.len());
        cascade.process_block_mono(chunk);
    }

    // Output must not be silent — the filter should ring
    let peak = buf.iter().map(|s| s.abs()).fold(0.0f32, f32::max);
    let energy: f32 = buf.iter().map(|s| s * s).sum();

    assert!(
        peak > 1e-6,
        "impulse response is silent (peak={peak}), filter coefficients may be broken"
    );
    assert!(
        energy > 1e-6,
        "impulse response has no energy ({energy}), filter is not processing"
    );

    // The filter should produce output beyond just the first sample
    // (a resonant filter rings for many samples)
    let tail_energy: f32 = buf[32..].iter().map(|s| s * s).sum();
    assert!(
        tail_energy > 1e-10,
        "filter has no tail ring (tail_energy={tail_energy}), not resonant"
    );
}

#[test]
fn forge_body_noise_produces_shaped_output() {
    let cart = load_acid_squelch();
    let mut cascade = Cascade::new();

    let coeffs = cart.interpolate(0.5, 0.5);

    // Feed white noise (deterministic seed via simple LCG)
    let mut buf = [0.0f32; 512];
    let mut seed: u32 = 12345;
    for s in buf.iter_mut() {
        seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
        *s = ((seed >> 16) as f32 / 32768.0) - 1.0;
    }

    let input_energy: f32 = buf.iter().map(|s| s * s).sum();
    for chunk in buf.chunks_mut(trench_core::BLOCK_SIZE) {
        cascade.set_targets(&coeffs, chunk.len());
        cascade.process_block_mono(chunk);
    }
    let output_energy: f32 = buf.iter().map(|s| s * s).sum();

    assert!(
        output_energy > 1e-6,
        "noise through forge body produced silence (energy={output_energy})"
    );

    // Output energy should differ from input — the filter is shaping the spectrum
    let ratio = output_energy / input_energy;
    assert!(
        (ratio - 1.0).abs() > 0.01,
        "output energy ratio={ratio:.4}, filter appears to be passthrough"
    );
}

#[test]
fn morph_changes_output() {
    let cart = load_acid_squelch();

    // Process same noise at morph=0 and morph=1, confirm different output
    let make_noise = || -> Vec<f32> {
        let mut buf = vec![0.0f32; 256];
        let mut seed: u32 = 99999;
        for s in buf.iter_mut() {
            seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
            *s = ((seed >> 16) as f32 / 32768.0) - 1.0;
        }
        buf
    };

    let mut buf_m0 = make_noise();
    let mut cascade_m0 = Cascade::new();
    let coeffs_m0 = cart.interpolate(0.0, 0.0);
    for chunk in buf_m0.chunks_mut(trench_core::BLOCK_SIZE) {
        cascade_m0.set_targets(&coeffs_m0, chunk.len());
        cascade_m0.process_block_mono(chunk);
    }

    let mut buf_m1 = make_noise();
    let mut cascade_m1 = Cascade::new();
    let coeffs_m1 = cart.interpolate(1.0, 0.0);
    for chunk in buf_m1.chunks_mut(trench_core::BLOCK_SIZE) {
        cascade_m1.set_targets(&coeffs_m1, chunk.len());
        cascade_m1.process_block_mono(chunk);
    }

    // The two outputs must be different
    let diff_energy: f32 = buf_m0
        .iter()
        .zip(buf_m1.iter())
        .map(|(a, b)| (a - b) * (a - b))
        .sum();

    assert!(
        diff_energy > 1e-4,
        "morph=0 and morph=1 produced identical output (diff_energy={diff_energy})"
    );
}
