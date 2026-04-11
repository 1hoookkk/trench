"""
TRENCH thumbwheel filmstrip renderer.
Run: blender --background --python tools/render_thumbwheel.py

E-mu style recessed horizontal thumbwheel on transparency.
Explicit parts only:
- outer drum halves (rotate)
- groove shoulder walls (rotate)
- inner emissive band (STATIONARY — light packet slides via UV offset)

No housing meshes.
No boolean cutters.
No topology-derived slit placement.
"""

import bmesh
import bpy
import math
import os
import shutil
import tempfile


STRIP_W = 123
STRIP_H = 28
FRAMES = 128

GLOW_CORE = (0.224, 0.878, 0.878, 1.0)   # #39E0E0
GLOW_BODY = (0.090, 0.749, 0.765, 1.0)   # #17BFC3
GLOW_TAIL = (0.047, 0.435, 0.451, 1.0)   # #0C6F73

CHARCOAL = (0.065, 0.058, 0.052, 1.0)    # lifted slightly for edge readability
GROOVE_BLACK = (0.025, 0.022, 0.020, 1.0)

OUT_DIR = os.path.join(
    os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
    "trench-plugin",
    "assets",
)

SEGMENTS = 128
RADIUS = 2.0
DRUM_LEN = 0.40
SHOULDER_LEN = 0.05         # reduced from 0.08
BAND_WIDTH = 0.19           # increased from 0.12
BAND_RADIUS = 1.88          # increased from 1.84 — slit reads more strongly
KNURL_RECESS = 0.030        # reduced from 0.045


def clear_scene():
    bpy.ops.wm.read_factory_settings(use_empty=True)
    scene = bpy.context.scene
    scene.render.film_transparent = True
    return scene


def make_drum_half(name, z_center):
    bpy.ops.mesh.primitive_cylinder_add(
        vertices=SEGMENTS,
        radius=RADIUS,
        depth=DRUM_LEN,
        location=(0.0, 0.0, z_center),
    )
    obj = bpy.context.active_object
    obj.name = name
    bpy.ops.object.shade_smooth()

    bpy.ops.object.mode_set(mode="EDIT")
    mesh = bmesh.from_edit_mesh(obj.data)

    for face in mesh.faces:
        face.select = False

    side_faces = [face for face in mesh.faces if len(face.verts) == 4]
    side_faces.sort(
        key=lambda face: math.atan2(
            face.calc_center_median().y,
            face.calc_center_median().x,
        )
    )

    channels = []
    for index, face in enumerate(side_faces):
        if index % 2 == 0:
            face.select = True
            channels.append(face)

    inset = (2.0 * math.pi * RADIUS / SEGMENTS) * 0.11
    bmesh.ops.inset_individual(mesh, faces=channels, thickness=inset)
    for face in channels:
        bmesh.ops.translate(mesh, vec=-face.normal * KNURL_RECESS, verts=face.verts)

    bmesh.update_edit_mesh(obj.data)
    bpy.ops.object.mode_set(mode="OBJECT")

    bevel = obj.modifiers.new("Bevel", "BEVEL")
    bevel.width = 0.012
    bevel.segments = 3
    bevel.profile = 0.45

    return obj


def make_shoulder(name, z_center, radius_low, radius_high):
    bpy.ops.mesh.primitive_cone_add(
        vertices=96,
        radius1=radius_low,
        radius2=radius_high,
        depth=SHOULDER_LEN,
        location=(0.0, 0.0, z_center),
    )
    obj = bpy.context.active_object
    obj.name = name
    bpy.ops.object.shade_smooth()
    return obj


def set_cylindrical_uv_z(obj, depth):
    mesh = obj.data
    uv_layer = mesh.uv_layers.new(name="UVMap")
    half_depth = depth * 0.5

    for polygon in mesh.polygons:
        for loop_index in polygon.loop_indices:
            loop = mesh.loops[loop_index]
            co = mesh.vertices[loop.vertex_index].co
            u = (math.atan2(co.y, co.x) / (2.0 * math.pi)) % 1.0
            v = 0.5 if half_depth == 0.0 else ((co.z / depth) + 0.5)
            uv_layer.data[loop_index].uv = (u, max(0.0, min(1.0, v)))


def build_materials():
    mat_drum = bpy.data.materials.new("DrumPlastic")
    mat_drum.use_nodes = True
    drum_nodes = mat_drum.node_tree.nodes
    drum_links = mat_drum.node_tree.links
    drum_nodes.clear()

    drum_out = drum_nodes.new("ShaderNodeOutputMaterial")
    drum_bsdf = drum_nodes.new("ShaderNodeBsdfPrincipled")
    drum_bsdf.inputs["Base Color"].default_value = CHARCOAL
    drum_bsdf.inputs["Metallic"].default_value = 0.0
    drum_bsdf.inputs["Roughness"].default_value = 0.54
    if "Specular IOR Level" in drum_bsdf.inputs:
        drum_bsdf.inputs["Specular IOR Level"].default_value = 0.28  # slightly more specular

    drum_noise = drum_nodes.new("ShaderNodeTexNoise")
    drum_noise.inputs["Scale"].default_value = 900.0
    drum_noise.inputs["Detail"].default_value = 12.0
    drum_noise.inputs["Roughness"].default_value = 0.65
    drum_bump = drum_nodes.new("ShaderNodeBump")
    drum_bump.inputs["Strength"].default_value = 0.12
    drum_bump.inputs["Distance"].default_value = 0.002

    drum_links.new(drum_noise.outputs["Fac"], drum_bump.inputs["Height"])
    drum_links.new(drum_bump.outputs["Normal"], drum_bsdf.inputs["Normal"])
    drum_links.new(drum_bsdf.outputs["BSDF"], drum_out.inputs["Surface"])

    mat_groove = bpy.data.materials.new("GrooveWall")
    mat_groove.use_nodes = True
    groove_nodes = mat_groove.node_tree.nodes
    groove_links = mat_groove.node_tree.links
    groove_nodes.clear()

    groove_out = groove_nodes.new("ShaderNodeOutputMaterial")
    groove_bsdf = groove_nodes.new("ShaderNodeBsdfPrincipled")
    groove_bsdf.inputs["Base Color"].default_value = GROOVE_BLACK
    groove_bsdf.inputs["Metallic"].default_value = 0.0
    groove_bsdf.inputs["Roughness"].default_value = 0.68
    if "Specular IOR Level" in groove_bsdf.inputs:
        groove_bsdf.inputs["Specular IOR Level"].default_value = 0.12
    groove_links.new(groove_bsdf.outputs["BSDF"], groove_out.inputs["Surface"])

    # ── Inner emissive band material ────────────────────────────────────
    mat_band = bpy.data.materials.new("InnerBand")
    mat_band.use_nodes = True
    band_nodes = mat_band.node_tree.nodes
    band_links = mat_band.node_tree.links
    band_nodes.clear()

    band_out = band_nodes.new("ShaderNodeOutputMaterial")
    band_bsdf = band_nodes.new("ShaderNodeBsdfPrincipled")
    band_bsdf.inputs["Base Color"].default_value = GROOVE_BLACK
    band_bsdf.inputs["Metallic"].default_value = 0.0
    band_bsdf.inputs["Roughness"].default_value = 0.62
    if "Specular IOR Level" in band_bsdf.inputs:
        band_bsdf.inputs["Specular IOR Level"].default_value = 0.12

    band_emit = band_nodes.new("ShaderNodeEmission")

    tex_coord = band_nodes.new("ShaderNodeTexCoord")
    separate_uv = band_nodes.new("ShaderNodeSeparateXYZ")
    band_links.new(tex_coord.outputs["UV"], separate_uv.inputs["Vector"])

    # Animated glow center — driven per-frame, band stays stationary
    glow_center = band_nodes.new("ShaderNodeValue")
    glow_center.name = "GlowCenter"
    glow_center.outputs[0].default_value = 0.25

    u_diff = band_nodes.new("ShaderNodeMath")
    u_diff.operation = "SUBTRACT"
    band_links.new(separate_uv.outputs["X"], u_diff.inputs[0])
    band_links.new(glow_center.outputs["Value"], u_diff.inputs[1])

    u_wrap = band_nodes.new("ShaderNodeMath")
    u_wrap.operation = "WRAP"
    u_wrap.inputs[1].default_value = -0.5
    u_wrap.inputs[2].default_value = 0.5
    band_links.new(u_diff.outputs["Value"], u_wrap.inputs[0])

    u_sign = band_nodes.new("ShaderNodeMath")
    u_sign.operation = "SIGN"
    band_links.new(u_wrap.outputs["Value"], u_sign.inputs[0])

    u_abs = band_nodes.new("ShaderNodeMath")
    u_abs.operation = "ABSOLUTE"
    band_links.new(u_wrap.outputs["Value"], u_abs.inputs[0])

    asymmetry = band_nodes.new("ShaderNodeMapRange")
    asymmetry.inputs[1].default_value = -1.0
    asymmetry.inputs[2].default_value = 1.0
    asymmetry.inputs[3].default_value = 0.72
    asymmetry.inputs[4].default_value = 1.35
    band_links.new(u_sign.outputs["Value"], asymmetry.inputs["Value"])

    u_shaped = band_nodes.new("ShaderNodeMath")
    u_shaped.operation = "MULTIPLY"
    band_links.new(u_abs.outputs["Value"], u_shaped.inputs[0])
    band_links.new(asymmetry.outputs["Result"], u_shaped.inputs[1])

    # Tighter, shorter glow packet
    core = band_nodes.new("ShaderNodeMapRange")
    core.interpolation_type = "SMOOTHSTEP"
    core.inputs[1].default_value = 0.0
    core.inputs[2].default_value = 0.020    # was 0.035 — tighter hot core
    core.inputs[3].default_value = 1.0
    core.inputs[4].default_value = 0.0
    band_links.new(u_shaped.outputs["Value"], core.inputs["Value"])

    body = band_nodes.new("ShaderNodeMapRange")
    body.interpolation_type = "SMOOTHSTEP"
    body.inputs[1].default_value = 0.0
    body.inputs[2].default_value = 0.045    # was 0.075 — shorter body
    body.inputs[3].default_value = 0.4
    body.inputs[4].default_value = 0.0
    band_links.new(u_shaped.outputs["Value"], body.inputs["Value"])

    tail = band_nodes.new("ShaderNodeMapRange")
    tail.interpolation_type = "SMOOTHSTEP"
    tail.inputs[1].default_value = 0.0
    tail.inputs[2].default_value = 0.08     # was 0.14 — shorter tail
    tail.inputs[3].default_value = 0.10
    tail.inputs[4].default_value = 0.0
    band_links.new(u_shaped.outputs["Value"], tail.inputs["Value"])

    glow_sum_a = band_nodes.new("ShaderNodeMath")
    glow_sum_a.operation = "ADD"
    band_links.new(core.outputs["Result"], glow_sum_a.inputs[0])
    band_links.new(body.outputs["Result"], glow_sum_a.inputs[1])

    glow_sum_b = band_nodes.new("ShaderNodeMath")
    glow_sum_b.operation = "ADD"
    band_links.new(glow_sum_a.outputs["Value"], glow_sum_b.inputs[0])
    band_links.new(tail.outputs["Result"], glow_sum_b.inputs[1])

    glow_clamp = band_nodes.new("ShaderNodeMath")
    glow_clamp.operation = "MINIMUM"
    glow_clamp.inputs[1].default_value = 1.0
    band_links.new(glow_sum_b.outputs["Value"], glow_clamp.inputs[0])

    upper_bias = band_nodes.new("ShaderNodeMapRange")
    upper_bias.interpolation_type = "SMOOTHSTEP"
    upper_bias.inputs[1].default_value = 0.0
    upper_bias.inputs[2].default_value = 1.0
    upper_bias.inputs[3].default_value = 0.42
    upper_bias.inputs[4].default_value = 1.0
    band_links.new(separate_uv.outputs["Y"], upper_bias.inputs["Value"])

    lower_bleed = band_nodes.new("ShaderNodeMapRange")
    lower_bleed.interpolation_type = "SMOOTHSTEP"
    lower_bleed.inputs[1].default_value = 0.0
    lower_bleed.inputs[2].default_value = 1.0
    lower_bleed.inputs[3].default_value = 0.16
    lower_bleed.inputs[4].default_value = 0.0
    band_links.new(separate_uv.outputs["Y"], lower_bleed.inputs["Value"])

    vertical_mix = band_nodes.new("ShaderNodeMath")
    vertical_mix.operation = "ADD"
    band_links.new(upper_bias.outputs["Result"], vertical_mix.inputs[0])
    band_links.new(lower_bleed.outputs["Result"], vertical_mix.inputs[1])

    vertical_clamp = band_nodes.new("ShaderNodeMath")
    vertical_clamp.operation = "MINIMUM"
    vertical_clamp.inputs[1].default_value = 1.0
    band_links.new(vertical_mix.outputs["Value"], vertical_clamp.inputs[0])

    final_glow = band_nodes.new("ShaderNodeMath")
    final_glow.operation = "MULTIPLY"
    band_links.new(glow_clamp.outputs["Value"], final_glow.inputs[0])
    band_links.new(vertical_clamp.outputs["Value"], final_glow.inputs[1])

    emit_strength = band_nodes.new("ShaderNodeMath")
    emit_strength.operation = "MULTIPLY"
    emit_strength.inputs[1].default_value = 30.0
    band_links.new(final_glow.outputs["Value"], emit_strength.inputs[0])
    band_links.new(emit_strength.outputs["Value"], band_emit.inputs["Strength"])

    color_ramp = band_nodes.new("ShaderNodeValToRGB")
    color_ramp.color_ramp.elements[0].position = 0.0
    color_ramp.color_ramp.elements[0].color = GLOW_CORE
    color_ramp.color_ramp.elements[1].position = 0.45
    color_ramp.color_ramp.elements[1].color = GLOW_BODY
    tail_stop = color_ramp.color_ramp.elements.new(1.0)
    tail_stop.color = GLOW_TAIL

    dist_norm = band_nodes.new("ShaderNodeMapRange")
    dist_norm.inputs[1].default_value = 0.0
    dist_norm.inputs[2].default_value = 0.08    # was 0.14 — matches shorter packet
    dist_norm.inputs[3].default_value = 0.0
    dist_norm.inputs[4].default_value = 1.0
    band_links.new(u_shaped.outputs["Value"], dist_norm.inputs["Value"])
    band_links.new(dist_norm.outputs["Result"], color_ramp.inputs["Fac"])
    band_links.new(color_ramp.outputs["Color"], band_emit.inputs["Color"])

    mix_shader = band_nodes.new("ShaderNodeMixShader")
    band_links.new(final_glow.outputs["Value"], mix_shader.inputs["Fac"])
    band_links.new(band_bsdf.outputs["BSDF"], mix_shader.inputs[1])
    band_links.new(band_emit.outputs["Emission"], mix_shader.inputs[2])
    band_links.new(mix_shader.outputs["Shader"], band_out.inputs["Surface"])

    return mat_drum, mat_groove, mat_band, glow_center


def configure_camera(scene):
    cam_data = bpy.data.cameras.new("Camera")
    cam_data.type = "ORTHO"
    cam_data.ortho_scale = 4.3      # tightened from 5.15
    cam_data.sensor_fit = "HORIZONTAL"

    cam = bpy.data.objects.new("Camera", cam_data)
    scene.collection.objects.link(cam)
    scene.camera = cam
    cam.location = (0.0, -7.0, 0.0)
    cam.rotation_euler = (math.radians(90.0), 0.0, 0.0)
    return cam


def configure_lights(scene):
    fill_data = bpy.data.lights.new("Fill", "AREA")
    fill_data.shape = "RECTANGLE"
    fill_data.size = 5.5
    fill_data.size_y = 0.7
    fill_data.energy = 14.0     # lifted from 11 for edge readability
    fill_data.color = (0.95, 0.92, 0.88)
    fill = bpy.data.objects.new("Fill", fill_data)
    scene.collection.objects.link(fill)
    fill.location = (0.0, -3.8, 1.2)
    fill.rotation_euler = (math.radians(74.0), 0.0, 0.0)

    rim_data = bpy.data.lights.new("Rim", "AREA")
    rim_data.shape = "RECTANGLE"
    rim_data.size = 4.0
    rim_data.size_y = 0.5
    rim_data.energy = 3.0       # lifted from 2.2
    rim_data.color = (0.82, 0.86, 0.88)
    rim = bpy.data.objects.new("Rim", rim_data)
    scene.collection.objects.link(rim)
    rim.location = (0.0, -3.0, -1.0)
    rim.rotation_euler = (math.radians(106.0), 0.0, 0.0)


def configure_render(scene):
    scene.render.engine = "CYCLES"
    scene.cycles.device = "GPU"
    scene.cycles.samples = 96
    scene.cycles.use_denoising = True
    scene.render.resolution_x = STRIP_W
    scene.render.resolution_y = STRIP_H
    scene.render.resolution_percentage = 100
    scene.render.image_settings.file_format = "PNG"
    scene.render.image_settings.color_mode = "RGBA"
    scene.render.film_transparent = True
    scene.frame_start = 1
    scene.frame_end = FRAMES


scene = clear_scene()

band_depth = BAND_WIDTH * 0.96
upper_drum_z = (BAND_WIDTH * 0.5) + SHOULDER_LEN + (DRUM_LEN * 0.5)
lower_drum_z = -upper_drum_z
upper_shoulder_z = (BAND_WIDTH * 0.5) + (SHOULDER_LEN * 0.5)
lower_shoulder_z = -upper_shoulder_z

drum_top = make_drum_half("DrumTop", upper_drum_z)
drum_bottom = make_drum_half("DrumBottom", lower_drum_z)
groove_upper = make_shoulder("GrooveUpper", upper_shoulder_z, BAND_RADIUS, RADIUS)
groove_lower = make_shoulder("GrooveLower", lower_shoulder_z, RADIUS, BAND_RADIUS)

bpy.ops.mesh.primitive_cylinder_add(
    vertices=96,
    radius=BAND_RADIUS,
    depth=band_depth,
    location=(0.0, 0.0, 0.0),
)
emissive_band = bpy.context.active_object
emissive_band.name = "EmissiveBand"
bpy.ops.object.shade_smooth()
set_cylindrical_uv_z(emissive_band, band_depth)

mat_drum, mat_groove, mat_band, glow_center = build_materials()
drum_top.data.materials.append(mat_drum)
drum_bottom.data.materials.append(mat_drum)
groove_upper.data.materials.append(mat_groove)
groove_lower.data.materials.append(mat_groove)
emissive_band.data.materials.append(mat_band)

# ── Animation ───────────────────────────────────────────────────────────
# ONLY the mechanical parts rotate. Emissive band stays stationary.
tooth_pitch = (2.0 * math.pi) / SEGMENTS
total_travel = math.pi * 0.20       # reduced from 0.35 — less rotation at this scale
rotation_rate = total_travel / (FRAMES - 1)
rotation_start = -0.5 * total_travel

for obj in [drum_top, drum_bottom, groove_upper, groove_lower]:
    obj.rotation_mode = "XYZ"
    driver = obj.driver_add("rotation_euler", 2)
    driver.driver.expression = f"{rotation_start} + {rotation_rate} * (frame - 1)"

# Emissive band does NOT rotate.
# Light packet slides via glow_center driven across UV space.
# Map frame 1→128 to UV U range that sweeps visible arc.
# The visible arc at the front of the cylinder is roughly U ∈ [0.25, 0.75].
# Sweep the glow center across that range.
glow_start = 0.30
glow_end = 0.70
glow_rate = (glow_end - glow_start) / (FRAMES - 1)

fc_glow = glow_center.outputs[0].driver_add("default_value")
fc_glow.driver.expression = f"{glow_start} + {glow_rate} * (frame - 1)"

configure_camera(scene)
configure_lights(scene)
configure_render(scene)

tmp = tempfile.mkdtemp(prefix="trench_wheel_")
print(f"Rendering {FRAMES} frames at {STRIP_W}x{STRIP_H} to {tmp}")

for frame in range(1, FRAMES + 1):
    scene.frame_set(frame)
    scene.render.filepath = os.path.join(tmp, f"frame_{frame:04d}.png")
    bpy.ops.render.render(write_still=True)
    if frame % 32 == 0:
        print(f"  {frame}/{FRAMES}")

print("Stitching filmstrip...")
strip_w = STRIP_W * FRAMES
strip = bpy.data.images.new("strip", strip_w, STRIP_H, alpha=True)
strip_px = [0.0] * (strip_w * STRIP_H * 4)

for index in range(FRAMES):
    frame_path = os.path.join(tmp, f"frame_{index + 1:04d}.png")
    image = bpy.data.images.load(frame_path)
    pixels = list(image.pixels)
    x_off = index * STRIP_W
    for y in range(STRIP_H):
        src = y * STRIP_W * 4
        dst = (y * strip_w + x_off) * 4
        strip_px[dst:dst + STRIP_W * 4] = pixels[src:src + STRIP_W * 4]
    bpy.data.images.remove(image)

strip.pixels = strip_px
os.makedirs(OUT_DIR, exist_ok=True)
out_path = os.path.join(OUT_DIR, "thumbwheel_strip.png")
strip.filepath_raw = out_path
strip.file_format = "PNG"
strip.save()
bpy.data.images.remove(strip)

shutil.rmtree(tmp, ignore_errors=True)
print(f"Done: {out_path}")
print(f"Filmstrip: {strip_w} x {STRIP_H} ({FRAMES} frames)")
