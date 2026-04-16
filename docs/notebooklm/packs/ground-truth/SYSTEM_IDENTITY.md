# TRENCH — System Identity

## What it is

TRENCH is a cartridge-based spectral weapon that gives producers shortcuts to taste through curated filter bodies, not open-ended design.

## What it is not

- Not a parametric EQ
- Not a filter design tool
- Not a workbench (that's TRENCHWORK, separate product)
- Not an emulation of any specific hardware

## Controls

Three. Morph, Q, cartridge selection. No drawing. No open-ended design surface ships with the plugin.

## The solver

The draw-the-filter solver is a **utility for internal authoring only**. It is not part of the product narrative. It is not visible to users. It uses RBJ cookbook formulas because curve-fitting is its job — it is not a character generation path. The solver does not ship as a user feature.

## Body taxonomy

| Lane | Description | Source |
|------|-------------|--------|
| Heritage Honest | Behaviors derived from heritage vocabulary, zeros preserved from donor material | P2K skins, Morpheus cubes via passthrough |
| Surgical Hybrid | Body-level zero-forcing applied to heritage pole skeletons | zero_law::apply on cube seeds |
| Betrayal | Coefficient-domain exploitation — stages that lie about what they are | Direct kernel-form authoring, morph-interior traps |

## Kill line

If a body can be replicated with a parametric EQ, it does not ship. The value is spectral character that requires 12-stage cascade interaction, zero-forcing, and betrayal to exist.

## Launch target

12 cartridges. Each has a name, a behavioral identity, a demo use case, and a reason it cannot be built with normal tools.

## Architecture (frozen)

- 12-stage serial DF2T cascade
- Kernel-form [c0..c4] interpolation only
- Q first, then morph interpolation order
- 32-sample control blocks, per-sample ramping
- Forge: 39062.5 Hz. Plugin: 44100 Hz
- Compiled-v1 (4-corner) cartridge format for plugin
- Stability by construction (r < 1), no pole sanitization
