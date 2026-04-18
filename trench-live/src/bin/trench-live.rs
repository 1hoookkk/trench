//! TRENCH LIVE standalone binary.
//!
//! Spawns the plugin with CPAL audio IO so the editor can be launched
//! outside a host for iteration and spike testing. The plugin itself
//! is the same code path as the VST3 / CLAP bundle — the standalone
//! wrapper is purely a host shim.
use nih_plug::prelude::*;
use trench_live::TrenchLive;

fn main() {
    nih_export_standalone::<TrenchLive>();
}
