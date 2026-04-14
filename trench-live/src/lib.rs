//! TRENCH LIVE — minimal Talking Hedz audition plugin.
//!
//! **Doctrine: the audio thread does nothing but math.**
//!
//! - No filesystem polling, no JSON reloading, no `Mutex`, no `RwLock`,
//!   no `Arc<AtomicCell<...>>` handoffs into `process()`.
//! - No editor — see `tools/bake_hedz_const.py` and
//!   `trench-core/src/hedz_rom.rs` for the provenance chain that
//!   produces the hardcoded Talking Hedz cartridge at build time.
//! - The cartridge is fully const at the Rust level (`trench_core::
//!   hedz_rom::HEDZ_CORNERS` + `HEDZ_BOOSTS`). `Cartridge::hedz_rom()`
//!   makes one `String` allocation for the name at plugin init —
//!   once — and nothing else ever allocates.
//! - `process()` walks 32-sample control blocks through
//!   `Cascade::set_targets` + `Cascade::set_boost` +
//!   `Cascade::process_block_mono`. Anything else added to this
//!   function is a doctrine violation.
//!
//! Parameter shape is intentionally small: two knobs (morph, Q) on
//! `[0, 1]`. Parameter reads on the hot path are `FloatParam::value()`,
//! which is lock-free per nih-plug's atomic parameter storage.
//!
//! If you find yourself wanting to live-reload a different cartridge:
//! rebuild the plugin with a different `hedz_rom.rs`. That is not a
//! workflow limitation, it is the point.

use nih_plug::prelude::*;
use std::sync::Arc;

use trench_core::{Cartridge, Cascade, BLOCK_SIZE};

pub struct TrenchLive {
    params: Arc<TrenchLiveParams>,
    cascade_l: Cascade,
    cascade_r: Cascade,
    cartridge: Cartridge,
}

#[derive(Params)]
pub struct TrenchLiveParams {
    #[id = "morph"]
    pub morph: FloatParam,
    #[id = "q"]
    pub q: FloatParam,
}

impl Default for TrenchLive {
    fn default() -> Self {
        Self {
            params: Arc::new(TrenchLiveParams::default()),
            cascade_l: Cascade::new(),
            cascade_r: Cascade::new(),
            // One allocation for the name String. Never touches the
            // audio thread. The 4×6×5 float cartridge body is copied
            // from `trench_core::hedz_rom::HEDZ_CORNERS` which lives
            // in .rodata.
            cartridge: Cartridge::hedz_rom(),
        }
    }
}

impl Default for TrenchLiveParams {
    fn default() -> Self {
        Self {
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

    // No editor — the plugin ships a hardcoded cartridge and exposes
    // two host-automation parameters. That's the whole UI surface.

    fn initialize(
        &mut self,
        _: &AudioIOLayout,
        buf: &BufferConfig,
        _: &mut impl InitContext<Self>,
    ) -> bool {
        // Hedz is a 44.1 kHz / 48 kHz body (see the heritage compile
        // path in `tools/bake_hedz_const.py`). Reject unusual rates
        // rather than silently alias the formants.
        let sr = buf.sample_rate;
        (sr - 44100.0).abs() <= 0.5 || (sr - 48000.0).abs() <= 0.5
    }

    fn reset(&mut self) {
        self.cascade_l.reset();
        self.cascade_r.reset();
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
        //    - call into trench_core::Cascade (pure DSP)
        //    - return ProcessStatus
        //  It must NOT:
        //    - take any lock of any kind
        //    - touch the filesystem
        //    - allocate (except via Cascade's stack-owned state)
        //    - parse JSON or any serialized format
        //    - block on channels, signals, or atomics that can spin
        //
        //  If you need to edit this function, re-read
        //  `pyruntime/CLAUDE.md` and `trench-core/CLAUDE.md` first.
        // ================================================================

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

            let mut channels = block.into_iter();
            if let Some(l) = channels.next() {
                self.cascade_l.process_block_mono(l);
            }
            if let Some(r) = channels.next() {
                self.cascade_r.process_block_mono(r);
            }

            if self.cascade_l.take_instability_flag() || self.cascade_r.take_instability_flag() {
                self.cascade_l.reset();
                self.cascade_r.reset();
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
