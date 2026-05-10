"""trap_vowel_slam - open dark 'uh' bracing to bright forward 'ih' slam.

10-number hand-authored body demonstrating the tools.body_dsl surface.
Lo: dark, mouth-back vowel; Hi: snapped-forward formant bite.
"""
from tools.body_dsl import Body, Slot

BODY = Body(
    name="trap_vowel_slam",
    boost=2.40,
    notes="dark uh braces into bright ih; mass+carve slam up an octave, pressure Q-climbs",
    slots=(
        Slot("anchor",   pole=(190, 170),   Q=0.94,            gain=(0.42, 0.40)),
        Slot("mass",     pole=(620, 1040),  Q=(0.980, 0.994),  gain=(0.58, 0.54)),
        Slot("carve",    pole=(1180, 2240), Q=(0.982, 0.996),  gain=(0.59, 0.55)),
        Slot("pressure", pole=(2280, 3040), Q=(0.992, 0.998),  gain=(0.60, 0.57)),
        Slot("detail",   pole=(4480, 5120), Q=(0.952, 0.980),  gain=(0.57, 0.52)),
        Slot("boundary", pole=(9400, 8420), Q=(0.974, 0.962),  gain=(0.52, 0.48)),
    ),
)
