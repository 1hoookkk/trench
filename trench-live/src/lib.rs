//! TRENCH LIVE — minimal Talking Hedz audition plugin.
//!
//! **Doctrine: the audio thread does nothing but math.**
//!
//! - No filesystem polling, no JSON reloading, no `Mutex`, no `RwLock`,
//!   no `Arc<AtomicCell<...>>` handoffs into `process()`.
//! - The vizia editor (`editor.rs`) runs on the GUI thread only and
//!   communicates with the audio thread exclusively through nih-plug's
//!   atomic parameter storage. Nothing in `editor.rs` may ever be
//!   invoked from `process()`.
//! - The shipping cartridge is hardcoded — see `tools/bake_hedz_const
//!   .py` and `trench-core/src/hedz_rom.rs` for the provenance chain
//!   that produces it at build time.
//! - `initialize()` calls `FilterEngine::prepare()` and loads the baked
//!   cartridge (one allocation for the name String, then nothing).
//! - `process()` calls `FilterEngine::process_block` on each channel
//!   slice. The engine handles down-resampling to NATIVE_SR (39062.5 Hz),
//!   the DF2T cascade, AGC, DC blocker, and up-resampling back to host SR.
//!   All scratch buffers are pre-allocated in `prepare()`; `process()` is
//!   allocation-free.
//! - The engine accepts any host sample rate. `initialize()` always returns
//!   `true`; resampling handles the conversion transparently.

use nih_plug::prelude::*;
use nih_plug_vizia::ViziaState;
use std::sync::Arc;

use trench_core::qsound_spatial::QSoundSpatial;
use trench_core::{Cartridge, FilterEngine};

mod editor;

pub struct TrenchLive {
    params: Arc<TrenchLiveParams>,
    engine_l: FilterEngine,
    engine_r: FilterEngine,
    spatial: QSoundSpatial,
    /// Host sample rate set in `initialize()`. Stored so `reset()` can
    /// re-call `reset_dsp()` without re-running the full `prepare()`.
    sample_rate: f64,
}

const HEDZ_SOURCE_JSON: &str =
    include_str!("../../trench-juce/plugin/assets/cartridges/P2k_013.json");

#[derive(Params)]
pub struct TrenchLiveParams {
    /// Persisted editor window state (size + scale factor). The GUI
    /// thread reads and writes this; the audio thread never touches it.
    #[persist = "editor-state"]
    pub editor_state: Arc<ViziaState>,

    #[id = "morph"]
    pub morph: FloatParam,
    #[id = "q"]
    pub q: FloatParam,
}

impl Default for TrenchLive {
    fn default() -> Self {
        Self {
            params: Arc::new(TrenchLiveParams::default()),
            engine_l: FilterEngine::new(),
            engine_r: FilterEngine::new(),
            spatial: QSoundSpatial::new(44_100.0),
            sample_rate: 44100.0,
        }
    }
}

impl Default for TrenchLiveParams {
    fn default() -> Self {
        Self {
            editor_state: editor::default_state(),

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

    fn params(&self) -> Arc<dyn Params> {
        self.params.clone()
    }

    fn editor(&mut self, _async_executor: AsyncExecutor<Self>) -> Option<Box<dyn Editor>> {
        editor::create(self.params.clone(), self.params.editor_state.clone())
    }

    fn initialize(
        &mut self,
        _: &AudioIOLayout,
        buf: &BufferConfig,
        _: &mut impl InitContext<Self>,
    ) -> bool {
        let sr = buf.sample_rate as f64;
        self.sample_rate = sr;
        // prepare() allocates scratch buffers and sinc tables — not on audio thread.
        self.engine_l.prepare(sr);
        self.engine_r.prepare(sr);
        self.spatial = QSoundSpatial::new(sr as f32);

        let cart = Cartridge::from_json(HEDZ_SOURCE_JSON).expect("embedded P2k_013 parses");
        if let Some(profile) = cart.spatial_profile.as_ref() {
            self.spatial.set_profile(profile);
            self.spatial.set_space(1.0);
        } else {
            self.spatial.set_space(0.0);
        }
        self.engine_l.load_cartridge(cart.clone());
        self.engine_r.load_cartridge(cart);
        true
    }

    fn reset(&mut self) {
        self.engine_l.reset_dsp();
        self.engine_r.reset_dsp();
        self.spatial.reset();
    }

    fn process(
        &mut self,
        buffer: &mut Buffer,
        _: &mut AuxiliaryBuffers,
        _: &mut impl ProcessContext<Self>,
    ) -> ProcessStatus {
        // ================================================================
        //  AUDIO-THREAD PURITY CHECKPOINT
        // ----------------------------------------------------------------
        //  This function must only:
        //    - read atomic parameter values (lock-free)
        //    - call into FilterEngine (pure DSP, allocation-free post-prepare)
        //    - return ProcessStatus
        //  FilterEngine::process_block uses pre-allocated scratch buffers
        //  and Option::take/put for borrow-checker safety — no heap calls.
        // ================================================================

        let morph = self.params.morph.value() as f64;
        let q = self.params.q.value() as f64;

        for (_, block) in buffer.iter_blocks(32) {
            let mut channels = block.into_iter();
            match (channels.next(), channels.next()) {
                (Some(l), Some(r)) => {
                    self.engine_l.process_block(l, morph, q);
                    self.engine_r.process_block(r, morph, q);
                    self.spatial.process_stereo(l, r);
                }
                (Some(l), None) => {
                    self.engine_l.process_block(l, morph, q);
                }
                (None, Some(r)) => {
                    self.engine_r.process_block(r, morph, q);
                }
                (None, None) => {}
            }

            if self.engine_l.take_instability_flag() || self.engine_r.take_instability_flag() {
                self.engine_l.reset_dsp();
                self.engine_r.reset_dsp();
                self.spatial.reset();
            }
        }

        ProcessStatus::Normal
    }
}

impl ClapPlugin for TrenchLive {
    const CLAP_ID: &'static str = "com.trenchwork.trench-live";
    const CLAP_DESCRIPTION: Option<&'static str> = Some("TRENCH Talking Hedz audition");
    const CLAP_MANUAL_URL: Option<&'static str> = None;
    const CLAP_SUPPORT_URL: Option<&'static str> = None;
    const CLAP_FEATURES: &'static [ClapFeature] = &[
        ClapFeature::AudioEffect,
        ClapFeature::Filter,
        ClapFeature::Stereo,
    ];
}

impl Vst3Plugin for TrenchLive {
    const VST3_CLASS_ID: [u8; 16] = *b"TrenchLiveAudPlu";
    const VST3_SUBCATEGORIES: &'static [Vst3SubCategory] =
        &[Vst3SubCategory::Fx, Vst3SubCategory::Filter];
}

nih_export_clap!(TrenchLive);
nih_export_vst3!(TrenchLive);
