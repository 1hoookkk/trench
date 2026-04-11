pub mod cartridge;
pub mod cascade;
pub mod cluster;
pub mod motor;
pub mod role;

pub use cartridge::{Cartridge, CornerData, NUM_COEFFS, NUM_STAGES};
pub use cascade::{Cascade, BLOCK_SIZE, TOTAL_STAGES};
pub use cluster::{Cluster, ClusterPlan, ClusterStage, InterferenceLink, StageFunction};
pub use motor::{ForceCurve, Motor, MotorId};
pub use role::Role;
