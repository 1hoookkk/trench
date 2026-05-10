# Peak/Shelf — research-only archive

Pre-veto authoring IRs and run logs from the previous Peak/Shelf authoring slice. Archived 2026-05-10.

**Do not use as a design surface.** TRENCH does not ship E-mu filters per `authoring/FILTER_WORK.md`. "Peak/Shelf" is the name of an E-mu Morpheus / EmulatorX3 filter class (DillusionMan tutorial: SHELF blends LP↔mid-shelf↔HP; PEAK is master volume; 06-order = 3 active stages). Authoring against this schema is a clone path and is shipping-blocked by FILTER_WORK.md vetos #1 and #4.

These files are kept only as:

- vocabulary reference for the SHELF / PEAK / FREQ controls,
- evidence of an authoring trajectory that was tried and rejected,
- input data for any future TRENCH-native body that wants to *measure distance* from Peak/Shelf as a reject-clones check.

Do NOT promote any file from this directory back into an active compiler, fixture, schema, or shipping cartridge path.

The compile path that produced these (`forge/compilers/compile_peak_shelf.py`, `peak_shelf_studio.py`, `schemas/peak_shelf_v1.json`, `runtime/trench-core/tests/peak_shelf_authoring.rs`, fixtures, audit/preview tools) was deleted in the same change that placed these files here. Cascade plot evidence at the time of deletion: the impl inverted the authored intent (M0 rolloff where intent asked for shelf lift; M100 narrow band-peak where intent asked for HF air). Both shipping-illegal and broken.

For the canonical heritage spec captured for reject-clones use, see the `project_peak_shelf_canonical_spec.md` memory.
