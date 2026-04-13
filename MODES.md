# TRENCH operating modes

## Shape Bank mode (DEFAULT)
- We collect shapes: static sonic locations + explicitly authored trajectories.
- The Shape Bank (`vault/_shapes`, 370+ bodies) is the Looperator token
  vocabulary.
- **Static targets are valid by definition.** Do not score, rank, prune, or
  "optimize them away".
- Job: correct artifacts and a truthful index, not curating taste.

## Trajectory mode (only when explicitly requested)
- When authoring a morph journey, judge the trajectory; a static plot is
  insufficient. Motion judgement applies only to trajectories.

## Operator mode (repo hygiene / reproducibility)
- Refactors, plumbing, indexing, fixing broken paths, reproducible runs.
- Prefer mechanical enforcement: single canonical truth files, regenerated
  inventories, stubs/mirrors.

## Shipping mode (release bodies)
- Follow `SPEC.md` + `BODIES.md`. Shipping bodies must pass the release gate.
- If Shape Bank rules conflict with Shipping rules, Shipping wins for shipping
  bodies.
