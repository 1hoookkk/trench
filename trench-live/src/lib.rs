use nih_plug::prelude::*;
use std::sync::Arc;
use std::time::Instant;

use crossbeam::atomic::AtomicCell;
use trench_core::{Cartridge, Cascade, BLOCK_SIZE};

mod editor;

fn live_json_path() -> std::path::PathBuf {
    let candidates = [
        std::env::var("TRENCH_LIVE_PATH").ok(),
        std::env::var("USERPROFILE")
            .ok()
            .map(|h| format!("{}/Trench/trench_live.json", h)),
    ];
    if let Some(c) = candidates.into_iter().flatten().next() {
        std::path::PathBuf::from(c)
    } else {
        std::path::PathBuf::from("trench_live.json")
    }
}

fn try_load_body(path: &std::path::Path) -> Option<Cartridge> {
    let json = std::fs::read_to_string(path).ok()?;
    Cartridge::from_json(&json).ok()
}

fn passthrough_cartridge() -> Cartridge {
    let identity_corner = [[1.0, 0.0, 0.0, 0.0, 0.0]; 6];
    Cartridge {
        name: "BYPASS".to_string(),
        corners: [identity_corner; 4],
        boosts: [1.0; 4],
    }
}

pub struct TrenchLive {
    params: Arc<TrenchLiveParams>,
    cascade_l: Cascade,
    cascade_r: Cascade,
    cartridge: Cartridge,
    cartridge_shared: Arc<std::sync::Mutex<Option<Cartridge>>>,
    body_name: Arc<AtomicCell<[u8; 64]>>,
    live_path: std::path::PathBuf,
    last_poll: Instant,
    last_modified: Option<std::time::SystemTime>,
    editor_state: Arc<editor::LiveEditorState>,
}

#[derive(Params)]
pub struct TrenchLiveParams {
    #[persist = "editor-state"]
    editor_state: Arc<editor::LiveEditorState>,

    #[id = "morph"]
    pub morph: FloatParam,
    #[id = "q"]
    pub q: FloatParam,
}

impl Default for TrenchLive {
    fn default() -> Self {
        let path = live_json_path();
        let cartridge = try_load_body(&path).unwrap_or_else(passthrough_cartridge);
        let mtime = std::fs::metadata(&path).ok().and_then(|m| m.modified().ok());

        let mut name_buf = [0u8; 64];
        let nb = cartridge.name.as_bytes();
        name_buf[..nb.len().min(63)].copy_from_slice(&nb[..nb.len().min(63)]);

        let params = Arc::new(TrenchLiveParams::default());
        let editor_state = params.editor_state.clone();

        Self {
            params,
            cascade_l: Cascade::new(),
            cascade_r: Cascade::new(),
            cartridge_shared: Arc::new(std::sync::Mutex::new(Some(cartridge.clone()))),
            body_name: Arc::new(AtomicCell::new(name_buf)),
            cartridge,
            live_path: path,
            last_poll: Instant::now(),
            last_modified: mtime,
            editor_state,
        }
    }
}

impl Default for TrenchLiveParams {
    fn default() -> Self {
        Self {
            editor_state: editor::LiveEditorState::new(),
            morph: FloatParam::new("Morph", 0.0, FloatRange::Linear { min: 0.0, max: 1.0 })
                .with_unit("%")
                .with_value_to_string(formatters::v2s_f32_percentage(0))
                .with_string_to_value(formatters::s2v_f32_percentage()),
            q: FloatParam::new("Q", 0.0, FloatRange::Linear { min: 0.0, max: 1.0 })
                .with_unit("%")
                .with_value_to_string(formatters::v2s_f32_percentage(0))
                .with_string_to_value(formatters::s2v_f32_percentage()),
        }
    }
}

impl TrenchLive {
    fn poll_file(&mut self) {
        if self.last_poll.elapsed().as_millis() < 500 { return; }
        self.last_poll = Instant::now();
        let mtime = match std::fs::metadata(&self.live_path) {
            Ok(m) => m.modified().ok(), Err(_) => return,
        };
        if mtime != self.last_modified {
            if let Some(cart) = try_load_body(&self.live_path) {
                nih_log!("TRENCH-LIVE: loaded '{}'", cart.name);
                let mut nb = [0u8; 64];
                let b = cart.name.as_bytes();
                nb[..b.len().min(63)].copy_from_slice(&b[..b.len().min(63)]);
                self.body_name.store(nb);
                if let Ok(mut g) = self.cartridge_shared.lock() { *g = Some(cart.clone()); }
                self.cartridge = cart;
                self.cascade_l.reset();
                self.cascade_r.reset();
            }
            self.last_modified = mtime;
        }
    }
}

impl Plugin for TrenchLive {
    const NAME: &'static str = "TRENCH LIVE";
    const VENDOR: &'static str = "Trenchwork";
    const URL: &'static str = "";
    const EMAIL: &'static str = "";
    const VERSION: &'static str = "0.1.0";
    const AUDIO_IO_LAYOUTS: &'static [AudioIOLayout] = &[AudioIOLayout {
        main_input_channels: NonZeroU32::new(2),
        main_output_channels: NonZeroU32::new(2),
        ..AudioIOLayout::const_default()
    }];
    const MIDI_INPUT: MidiConfig = MidiConfig::None;
    const SAMPLE_ACCURATE_AUTOMATION: bool = true;
    type SysExMessage = ();
    type BackgroundTask = ();

    fn params(&self) -> Arc<dyn Params> { self.params.clone() }

    fn editor(&mut self, _: AsyncExecutor<Self>) -> Option<Box<dyn Editor>> {
        Some(Box::new(editor::LiveEditor::new(
            self.params.clone(),
            self.editor_state.clone(),
            self.cartridge_shared.clone(),
            self.body_name.clone(),
        )))
    }

    fn initialize(&mut self, _: &AudioIOLayout, buf: &BufferConfig, _: &mut impl InitContext<Self>) -> bool {
        let sr = buf.sample_rate;
        (sr - 44100.0).abs() <= 0.5 || (sr - 48000.0).abs() <= 0.5
    }

    fn reset(&mut self) {
        self.cascade_l.reset();
        self.cascade_r.reset();
    }

    fn process(&mut self, buffer: &mut Buffer, _: &mut AuxiliaryBuffers, _: &mut impl ProcessContext<Self>) -> ProcessStatus {
        self.poll_file();
        let morph = self.params.morph.value() as f64;
        let q = self.params.q.value() as f64;
        let coeffs = self.cartridge.interpolate(morph, q);
        let boost = self.cartridge.interpolate_boost(morph, q);

        for (_, block) in buffer.iter_blocks(BLOCK_SIZE) {
            let n = block.samples();
            self.cascade_l.set_targets(&coeffs, n);
            self.cascade_r.set_targets(&coeffs, n);
            self.cascade_l.set_boost(boost, n);
            self.cascade_r.set_boost(boost, n);
            let mut ch = block.into_iter();
            if let Some(l) = ch.next() { self.cascade_l.process_block_mono(l); }
            if let Some(r) = ch.next() { self.cascade_r.process_block_mono(r); }
            if self.cascade_l.take_instability_flag() || self.cascade_r.take_instability_flag() {
                self.cascade_l.reset(); self.cascade_r.reset();
            }
        }
        ProcessStatus::Normal
    }
}

impl ClapPlugin for TrenchLive {
    const CLAP_ID: &'static str = "com.trenchwork.trench-live";
    const CLAP_DESCRIPTION: Option<&'static str> = Some("TRENCH live body audition");
    const CLAP_MANUAL_URL: Option<&'static str> = None;
    const CLAP_SUPPORT_URL: Option<&'static str> = None;
    const CLAP_FEATURES: &'static [ClapFeature] = &[ClapFeature::AudioEffect, ClapFeature::Filter, ClapFeature::Stereo];
}

impl Vst3Plugin for TrenchLive {
    const VST3_CLASS_ID: [u8; 16] = *b"TrenchLiveAudPlu";
    const VST3_SUBCATEGORIES: &'static [Vst3SubCategory] = &[Vst3SubCategory::Fx, Vst3SubCategory::Filter];
}

nih_export_clap!(TrenchLive);
nih_export_vst3!(TrenchLive);
