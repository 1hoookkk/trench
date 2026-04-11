# Better Than E-MU Holdout Scorecard

- Design corpus: `C:\Users\hooki\Trench\datasets\wav_corpus\design_split_emotion_strat`
- Holdout corpus: `C:\Users\hooki\Trench\datasets\wav_corpus\holdout_split_emotion_strat`
- Design probe mode: `corpus`
- Holdout probe mode: `corpus`
- Allow fallback: `False`
- Max abs metric contribution: `1.5`
- Balanced sampling: `True`
- Balance seed: `42`
- Max clips per class: `None`
- Design balance input counts: `{'bass': 56, 'vocal': 56, 'drum': 56, 'mix': 56}`
- Design balance output counts: `{'bass': 56, 'vocal': 56, 'drum': 56, 'mix': 56}`
- Design balance target count: `56`
- Holdout balance input counts: `{'bass': 14, 'vocal': 14, 'drum': 14, 'mix': 14}`
- Holdout balance output counts: `{'bass': 14, 'vocal': 14, 'drum': 14, 'mix': 14}`
- Holdout balance target count: `14`
- Design label counts: `{'bass': 56, 'vocal': 56, 'drum': 56, 'mix': 56}`
- Holdout label counts: `{'bass': 14, 'vocal': 14, 'drum': 14, 'mix': 14}`
- Design missing labels: `[]`
- Holdout missing labels: `[]`

## Summary

| Body | Design Winner | Holdout Promoted | Holdout Gate OK | Holdout Wins | Holdout Composite | Holdout Oracle |
|---|---|---:|---:|---:|---:|---|
| Speaker Knockerz | Speaker Knockerz C4 | True | True | 3 | 2.810 | Speaker Knockerz C4 |
| Aluminum Siding | Aluminum Siding C2 | True | True | 3 | 3.000 | Aluminum Siding C2 |
| Small Talk Ah-Ee | Small Talk Ah-Ee C2 | True | True | 4 | 3.996 | Small Talk Ah-Ee C2 |
| Cul-De-Sac | Cul-De-Sac C3 | False | False | 3 | 3.231 | Cul-De-Sac C3 |

## Rule

- Candidate selection is based on design corpus only.
- Final promotion decision is holdout-only.