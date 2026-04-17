//! TRENCH FORGE DEV — hot-reload audition harness for body authoring.
//!
//! # Scoped doctrine exception (dev-only, non-shipping)
//!
//! The shipping plugin (`trench-live`) bans filesystem polling, JSON
//! reloading, `Mutex`, `RwLock`, and `Arc<T>` handoffs on the audio
//! thread. **Those bans remain intact in shipping.** This crate carves
//! a narrowly scoped exception for the write-hear authoring loop:
//!
//! - **Non-shipping.** `trench-forge-dev` must never be promoted into a
//!   release build. Shipping is `trench-live` (fixed cartridge) and
//!   `trench-juce/plugin/` (JUCE standalone + VST3). This crate exists
//!   only to close the latency between `compile_raw.py` output and
//!   DAW audition during body authoring.
//! - **Audio thread is still pure.** `process()` only reads atomic
//!   parameter values, performs a wait-free `ArcSwap::load()` to pin
//!   the active `Cartridge`, and calls into `trench_core::Cascade`.
//!   No filesystem touch, no JSON parse, no allocation, no blocking
//!   of any kind.
//! - **Control thread handles IO.** A background `notify` watcher
//!   thread monitors the JSON file at `TRENCH_FORGE_WATCH` (default:
//!   `./forge-watch/active.json`). On change: read, parse via
//!   `Cartridge::from_json`, and `store()` into the shared
//!   `ArcSwap<Cartridge>`. On parse failure the previous cartridge
//!   is retained — malformed writes never leak into the DSP path.
//! - **Atomic handoff.** Coefficient transitions between cartridges
//!   are absorbed by the existing `Cascade::set_targets` ramp
//!   machinery — the same mechanism shipping uses for MORPH
//!   automation. No clicks, no instability.
//!
//! If the watch file is missing or invalid at startup, the harness
//! falls back to the const `Cartridge::hedz_rom()` so the plugin is
//! always auditable.

use std::num::NonZeroU32;
use std::path::PathBuf;
use std::sync::Arc;

use arc_swap::ArcSwap;
use nih_plug::prelude::*;
use notify::{EventKind, RecommendedWatcher, RecursiveMode, Watcher};

use trench_core::{Cartridge, Cascade, BLOCK_SIZE};

const DEFAULT_WATCH_PATH: &str = "./forge-watch/active.json";

pub struct TrenchForgeDev {
    params: Arc<TrenchForgeDevParams>,
    cascade_l: Cascade,
    cascade_r: Cascade,
    cartridge: Arc<ArcSwap<Cartridge>>,
    // Kept alive so the watcher thread keeps running. Never touched by
    // the audio thread. Option<> because notify init may fail on
    // exotic filesystems; plugin still works as a fixed-cartridge
    // fallback in that case.
    _watcher: Option<RecommendedWatcher>,
}

#[derive(Params)]
pub struct TrenchForgeDevParams {
    #[id = "morph"]
    pub morph: FloatParam,
    #[id = "q"]
    pub q: FloatParam,
}

impl Default for TrenchForgeDev {
    fn default() -> Self {
        Self {
            params: Arc::new(TrenchForgeDevParams::default()),
            cascade_l: Cascade::new(),
            cascade_r: Cascade::new(),
            // Fallback cartridge so the plugin is always auditable
            // even before the first hot-reload lands.
            cartridge: Arc::new(ArcSwap::from(Arc::new(Cartridge::hedz_rom()))),
            _watcher: None,
        }
    }
}

impl Default for TrenchForgeDevParams {
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

fn watch_path() -> PathBuf {
    std::env::var("TRENCH_FORGE_WATCH")
        .map(PathBuf::from)
        .unwrap_or_else(|_| PathBuf::from(DEFAULT_WATCH_PATH))
}

/// Try to read + parse the watch file and swap the active cartridge.
/// Silent on failure — the previous cartridge is retained, which is the
/// correct behavior for transient IO errors or partial writes.
fn try_reload(path: &PathBuf, swap: &Arc<ArcSwap<Cartridge>>) {
    if let Ok(json) = std::fs::read_to_string(path) {
        if let Ok(cart) = Cartridge::from_json(&json) {
            swap.store(Arc::new(cart));
        }
    }
}

impl Plugin for TrenchForgeDev {
    const NAME: &'static str = "TRENCH FORGE DEV";
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

    fn initialize(
        &mut self,
        _: &AudioIOLayout,
        buf: &BufferConfig,
        _: &mut impl InitContext<Self>,
    ) -> bool {
        // Sample-rate gate — same as trench-live. Hedz-lineage bodies
        // are 44.1/48 kHz; reject anything else rather than silently
        // alias the formants.
        let sr = buf.sample_rate;
        if !((sr - 44100.0).abs() <= 0.5 || (sr - 48000.0).abs() <= 0.5) {
            return false;
        }

        // Set up the FS watcher on a background thread. All file IO
        // and JSON parsing happens here, never in process().
        let path = watch_path();
        let swap = self.cartridge.clone();

        // Initial load if the file already exists — so if the user
        // launches with a JSON already in place we pick it up on start.
        try_reload(&path, &swap);

        if let Some(dir) = path.parent().filter(|p| !p.as_os_str().is_empty() && p.exists())
        {
            let watch_target = path.clone();
            let swap_for_cb = swap.clone();
            let watcher_result: notify::Result<RecommendedWatcher> =
                notify::recommended_watcher(move |res: notify::Result<notify::Event>| {
                    if let Ok(event) = res {
                        // React to create/modify events that touch our file.
                        let touches_target = event.paths.iter().any(|p| p == &watch_target);
                        let relevant = matches!(
                            event.kind,
                            EventKind::Modify(_) | EventKind::Create(_)
                        );
                        if touches_target && relevant {
                            try_reload(&watch_target, &swap_for_cb);
                        }
                    }
                });

            if let Ok(mut watcher) = watcher_result {
                let _ = watcher.watch(dir, RecursiveMode::NonRecursive);
                self._watcher = Some(watcher);
            }
        }

        true
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
        //  AUDIO-THREAD PURITY CHECKPOINT (dev-harness scope)
        // ----------------------------------------------------------------
        //  This function must only:
        //    - read atomic parameter values (lock-free)
        //    - perform a wait-free ArcSwap::load() of the current
        //      Cartridge (lock-free; pin via the returned Guard)
        //    - call into trench_core::Cascade (pure DSP)
        //    - return ProcessStatus
        //  It must NOT:
        //    - take any lock of any kind (Mutex, RwLock, etc.)
        //    - touch the filesystem
        //    - allocate (except via Cascade's stack-owned state)
        //    - parse JSON or any serialized format
        //    - block on channels, signals, or atomics that can spin
        //
        //  ArcSwap::load() is wait-free per `arc-swap`'s contract —
        //  it's a stronger guarantee than Mutex and is the ONLY
        //  concession from the shipping trench-live purity rules.
        //  Anything else added to this function violates the
        //  dev-harness exception.
        // ================================================================

        let morph = self.params.morph.value() as f64;
        let q = self.params.q.value() as f64;

        let cart_guard = self.cartridge.load();
        let cart: &Cartridge = cart_guard.as_ref();
        let coeffs = cart.interpolate(morph, q);
        let boost = cart.interpolate_boost(morph, q);

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

impl ClapPlugin for TrenchForgeDev {
    const CLAP_ID: &'static str = "com.trenchwork.trench-forge-dev";
    const CLAP_DESCRIPTION: Option<&'static str> =
        Some("TRENCH hot-reload authoring harness (dev, non-shipping)");
    const CLAP_MANUAL_URL: Option<&'static str> = None;
    const CLAP_SUPPORT_URL: Option<&'static str> = None;
    const CLAP_FEATURES: &'static [ClapFeature] = &[
        ClapFeature::AudioEffect,
        ClapFeature::Filter,
        ClapFeature::Stereo,
    ];
}

impl Vst3Plugin for TrenchForgeDev {
    const VST3_CLASS_ID: [u8; 16] = *b"TrenchForgeDev01";
    const VST3_SUBCATEGORIES: &'static [Vst3SubCategory] =
        &[Vst3SubCategory::Fx, Vst3SubCategory::Filter];
}

nih_export_clap!(TrenchForgeDev);
nih_export_vst3!(TrenchForgeDev);
