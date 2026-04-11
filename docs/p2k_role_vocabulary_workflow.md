# P2k Role Vocabulary Workflow

Use P2k as **training data for stage-role vocabulary**, not as an audio metric baseline.

## 1) Build vocabulary from P2k_000..032

```bash
python tools/build_p2k_role_vocabulary.py
```

Output:
- `datasets/role_vocab/p2k_role_vocabulary_v1.json`

## 2) Build shipping name targets from mapped P2k bodies

```bash
python tools/build_shipping_role_targets.py
```

Output:
- `datasets/role_vocab/shipping_role_targets_v1.json`

## 3) Rank shipping finalists by role vocabulary distance

```bash
python tools/promote_by_role_vocabulary.py
```

Output:
- `vault/_scorecards/role_vocab_scorecard.json`

## 4) Forge directly against role vocabulary objective

```bash
python tools/forge_pymoo.py speaker_knockerz --quality-mode role_vocab --pop 80 --gen 60
python tools/forge_pymoo.py aluminum_siding --quality-mode role_vocab --pop 80 --gen 60
python tools/forge_pymoo.py small_talk_ah_ee --quality-mode role_vocab --pop 80 --gen 60
python tools/forge_pymoo.py cul_de_sac --quality-mode role_vocab --pop 80 --gen 60
```

Quality modes:
- `role_vocab` (default): optimize to P2k-trained stage-role targets.
- `legacy`: old audio-shape quality objective.
- `hybrid`: role-vocab objective + light legacy tie-break.

## Notes
- This workflow uses only P2k corpus for vocabulary training.
- Shipping names remain the only normative label anchor for target selection.
