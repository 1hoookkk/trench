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

When `verdicts.csv` exists, the sampler automatically applies a weighted centroid
shift toward keeps and away from rejects. 25% of the batch is drawn from the
unshifted envelope (exploration). See `tools/workbench/bias.py` for the v1 contract.

Options:

```
--scope           Scope name (default: hf_resonance)
--n               Target admitted count (default: 20)
--seed            RNG seed for reproducibility
--out-dir         Output directory (default: _queue/)
--no-bias         Disable centroid shift even if verdicts.csv exists
--fresh           Wipe existing scope pills from _queue/ before generating
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

To re-audition pills previously logged as 'maybe':

```
python tools/workbench/verdict.py --resume --revisit-maybes
```

---

## Step 3 — Promote keepers

Pills with verdict `keep` or `maybe` are candidates for the shipping body roster.
To promote one, copy it from `_queue/` to `cartridges/pills/` and run `./check`.

---

## Step 4 — Re-sample with bias

Run the sampler again. It reads `verdicts.csv` automatically and shifts the scope
envelope toward your keeps, away from your rejects, then refills the queue:

```
python tools/workbench/sampler.py --n 20
```

The bias report prints per-corner center shifts and exploration fraction.
To start completely fresh (wipe queue, ignore previous verdicts in envelope):

```
python tools/workbench/sampler.py --fresh --no-bias --n 20
```

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
| `sampler.py` | Generate admitted pills into `_queue/` (with bias) |
| `verdict.py` | Audition loop + CSV logger |
| `bias.py` | Weighted centroid shift — v1 contract |
| `_queue/` | Persistent candidates (gitignored) |
| `verdicts.csv` | Persistent verdict log (gitignored) |
