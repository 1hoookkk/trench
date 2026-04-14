# Synthetic NotebookLM filter template excerpt
#
# Used by tests/test_extract_emu_filter_params.py to pin the parser to
# the bullet format documented in the task description.
#
# === REAL DATA (do not edit) ===
# `synthetic_hedz` stage 1 is the exact (type, low-freq, low-gain,
# high-freq, high-gain) tuple the user quoted verbatim from the real
# Talking Hedz `13-templates-filter-04.md` NotebookLM export. The test
# locks these five integers — if you change them, the parser truth
# anchor for stage 1 is gone and the test will fail.
#
# === FABRICATED SCAFFOLDING (clearly fake) ===
# Stages 2-6 of `synthetic_hedz` and the entire `synthetic_acwow`
# template are NOT from any E-mu source. They use obviously synthetic
# monotonic patterns so nobody can mistake them for canonical heritage
# data. Their only purpose is to exercise:
#   - 6-section grouping per template
#   - multi-template parsing
#   - the type-absolute=108 path documented in
#     docs/looperator/morph_designer_behavior.md:48
#   - missing-section zero-fill (synthetic_acwow only sets sections 1
#     and 3; sections 2,4,5,6 must default to zero)
#
# Once the real `10-13-templates-filter-*.md` exports land in this
# directory, delete this synthetic fixture and update the test to read
# the real bundle.

- synthetic_hedz/filter/type-absolute: 144
- synthetic_hedz/filter/frequency: 0
- synthetic_hedz/filter/gain: 0
- synthetic_hedz/filter/designer-section[1]/type: 3
- synthetic_hedz/filter/designer-section[1]/low-freq: 0
- synthetic_hedz/filter/designer-section[1]/low-gain: 127
- synthetic_hedz/filter/designer-section[1]/high-freq: 116
- synthetic_hedz/filter/designer-section[1]/high-gain: 87
- synthetic_hedz/filter/designer-section[2]/type: 1
- synthetic_hedz/filter/designer-section[2]/low-freq: 2
- synthetic_hedz/filter/designer-section[2]/low-gain: 20
- synthetic_hedz/filter/designer-section[2]/high-freq: 22
- synthetic_hedz/filter/designer-section[2]/high-gain: 220
- synthetic_hedz/filter/designer-section[3]/type: 1
- synthetic_hedz/filter/designer-section[3]/low-freq: 3
- synthetic_hedz/filter/designer-section[3]/low-gain: 30
- synthetic_hedz/filter/designer-section[3]/high-freq: 33
- synthetic_hedz/filter/designer-section[3]/high-gain: 230
- synthetic_hedz/filter/designer-section[4]/type: 1
- synthetic_hedz/filter/designer-section[4]/low-freq: 4
- synthetic_hedz/filter/designer-section[4]/low-gain: 40
- synthetic_hedz/filter/designer-section[4]/high-freq: 44
- synthetic_hedz/filter/designer-section[4]/high-gain: 240
- synthetic_hedz/filter/designer-section[5]/type: 1
- synthetic_hedz/filter/designer-section[5]/low-freq: 5
- synthetic_hedz/filter/designer-section[5]/low-gain: 50
- synthetic_hedz/filter/designer-section[5]/high-freq: 55
- synthetic_hedz/filter/designer-section[5]/high-gain: 250
- synthetic_hedz/filter/designer-section[6]/type: 1
- synthetic_hedz/filter/designer-section[6]/low-freq: 6
- synthetic_hedz/filter/designer-section[6]/low-gain: 60
- synthetic_hedz/filter/designer-section[6]/high-freq: 66
- synthetic_hedz/filter/designer-section[6]/high-gain: 260

- synthetic_acwow/filter/type-absolute: 108
- synthetic_acwow/filter/frequency: 0.5
- synthetic_acwow/filter/gain: 0.25
- synthetic_acwow/filter/designer-section[1]/type: 2
- synthetic_acwow/filter/designer-section[1]/low-freq: 11
- synthetic_acwow/filter/designer-section[1]/low-gain: 111
- synthetic_acwow/filter/designer-section[1]/high-freq: 11
- synthetic_acwow/filter/designer-section[1]/high-gain: 111
- synthetic_acwow/filter/designer-section[3]/type: 2
- synthetic_acwow/filter/designer-section[3]/low-freq: 33
- synthetic_acwow/filter/designer-section[3]/low-gain: 133
- synthetic_acwow/filter/designer-section[3]/high-freq: 33
- synthetic_acwow/filter/designer-section[3]/high-gain: 133
