# Source
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_cube_display_model.json` (113314 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\analysis\build_filter_cube_display_model.py` (8063 bytes)
- `C:\Users\hooki\trench_re_vault\scratch\analysis\filter_cube_display_report.md` (2985 bytes)

# What this data is
This dataset represents the structural identity of the morph filter templates extracted from the original E-mu hardware. Instead of relying on file names (which are alias-heavy and sometimes duplicated), the reverse-engineering process parsed each filter's six internal `designer-section` layers. This yields a deterministic `shape_id` based on the topology of `type`, `low_gain`, `low_freq`, `high_gain`, and `high_freq` across the layers.

# Why it's here
The Looperator time-grid will visualize each placed token as a filter cube preview. To correctly draw these preview cubes and group equivalent templates, the CLEAN repo needs this mapping to distinguish unique structural topologies from mere preset aliases. 

# Sanitized content

## Shape Dictionary Schema
The core extraction distills 79 templates down to 72 unique shapes. Below is the parsed schema that Looperator should expect:

```json
{
  "summary": {
    "template_count": 79,
    "unique_shape_count": 72,
    "deduplicated_template_count": 7,
    "type_absolute_groups": {
      "108": { "template_count": 70, "unique_shape_count": 69 },
      "144": { "template_count": 9, "unique_shape_count": 3 }
    }
  },
  "shapes": [
    {
      "shape_id": "string",
      "type_absolute": "integer",
      "display_name": "string",
      "alias_count": "integer",
      "aliases": [
        {
          "template_name": "string",
          "file": "string",
          "frequency_default": "float",
          "gain_default": "float"
        }
      ],
      "frequency_defaults": ["float"],
      "gain_defaults": ["float"],
      "designer_sections": [
        {
          "index": "integer (1 to 6)",
          "type": "integer",
          "low_gain": "integer (0-127)",
          "low_freq": "integer (0-127)",
          "high_gain": "integer (0-127)",
          "high_freq": "integer (0-127)"
        }
      ]
    }
  ]
}
```

## Collapsed Alias Groups
Certain filter presets are structurally identical. Looperator should present these as aliases under a single visual cube:
- **`abs108_b352f374685b`** (Type 108): `AC Wow`, `AC Morph`
- **`abs144_b5c825095c87`** (Type 144): `50`, `100`, `hedz0`, `d`, `hedz`
- **`abs144_7949e309f941`** (Type 144): `morph0q100`, `morph100q100`, `ssssss`

# Integration notes
- Store tokens in the CLEAN repo using `shape_id` as the primary invariant key.
- The `type_absolute` field should only be used as a secondary grouping layer (e.g., separating 108s from 144s in the picker).
- Display the `display_name` property as the main UI label. Treat the other names under `aliases` strictly as synonymous preset loader options.
- Maintain `frequency_defaults` and `gain_defaults` as preset parameters that a user can snap to, rather than treating them as part of the immutable response curve.
- The final rendering loop must interpret the six `designer_sections` and visualize those parameters, rather than reading hard-coded arrays.

# UNKNOWNs
- **Coordinate System & Geometry:** The RE extraction parses the six sections (`type`, `low/high gain`, `low/high freq`) but does not explicitly define how these five parameters map onto the 4-corner morph-×-Q coordinate space, nor does it define rendering bounds for thumbnails.
- **Color/Channel Encoding:** Missing mapping data on how section channels or types should be color-differentiated in the cube UI.
- **Adjacency Hints:** Missing any explicit proximity links. (Likely handled by the separate Atlas data, but worth noting in relation to display).
