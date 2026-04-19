# Workbench loop

Rapid body triage: generate pill candidates → audition in the live plugin → log verdicts.

---

## Prerequisites

- TRENCH plugin running in your DAW (or standalone), polling `TRENCH_LIVE_PATH`
- `pip install numpy` (pill_gate dependency)

Default live path (no env var needed):

```
~/Trench/trench_live.json
```

Set `TRENCH_LIVE_PATH` to override.

---

## Step 1 — Fill the queue

Generate 20 admitted pill candidates for the `hf_resonance` scope:

```
python tools/workbench/sampler.py --scope hf_resonance --n 20
```

Candidates land in `tools/workbench/_queue/`. Only gate-passing pills are written.
Rejected candidates print their fail reason and are discarded.

Options:

```
--scope     Scope name (default: hf_resonance)
--n         Target admitted count (default: 20)
--seed      RNG seed for reproducibility
--out-dir   Output directory (default: _queue/)
```

---

## Step 2 — Audition and log verdicts

```
python tools/workbench/verdict.py
```

Each pill is hot-swapped into the running plugin. Press one key per pill:

```
1 = reject    2 = maybe    3 = keep    q/Esc = quit
```

Verdicts append to `tools/workbench/verdicts.csv`.

To resume a previous session (skip already-logged pills):

```
python tools/workbench/verdict.py --resume
```

---

## Step 3 — Promote keepers

Pills with verdict `keep` or `maybe` are candidates for the shipping body roster.
To promote one, copy it from `_queue/` to `cartridges/pills/` and run `./check`.

---

## Step 4 — Re-sample (coming)

A bias re-sampler will read `verdicts.csv` and re-weight the scope envelope toward
high-rated candidates. Design checkpoint with user before implementing.

---

## Scope caveat

`SCOPE_ENVELOPES` in `tools/cube_authoring/reroll.py` has one entry: `hf_resonance`.
Each new body scope needs an envelope authored there before this loop can generate
candidates for it. Scope authoring is a product decision — do not add scopes without
user direction.

---

## Files

| Path | Purpose |
|------|---------|
| `sampler.py` | Generate admitted pills into `_queue/` |
| `verdict.py` | Audition loop + CSV logger |
| `_queue/` | Ephemeral candidates (gitignored) |
| `verdicts.csv` | Persistent verdict log (gitignored) |
