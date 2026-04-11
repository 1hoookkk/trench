"""Run in Blender with the thumbwheel scene already loaded.
Changes material to bone white plastic, sets render size to 135x35 (UI half-scale),
and renders 128 frames to trench-plugin/assets/thumbwheel/.

Usage: In Blender scripting tab, open and run this file.
Or: blender thumbwheel_setup.blend --background --python tools/patch_bone_white.py
"""
import bpy
import os

# --- 1. Patch material to bone white plastic ---
mat = None
for m in bpy.data.materials:
    if 'Copper' in m.name or 'Matte' in m.name or 'Wheel' in m.name:
        mat = m
        break

if mat is None:
    # Try the first material on the Thumbwheel object
    wheel = bpy.data.objects.get("Thumbwheel")
    if wheel and wheel.data.materials:
        mat = wheel.data.materials[0]

if mat is None:
    print("ERROR: No material found. Ensure Thumbwheel object exists with a material.")
else:
    print(f"Patching material: {mat.name}")
    nodes = mat.node_tree.nodes

    # Find the Principled BSDF
    principled = None
    for n in nodes:
        if n.type == 'BSDF_PRINCIPLED':
            principled = n
            break

    if principled:
        # Bone white plastic: warm off-white, non-metallic, slight roughness
        principled.inputs['Base Color'].default_value = (0.85, 0.82, 0.78, 1.0)  # warm bone
        principled.inputs['Metallic'].default_value = 0.0  # plastic, not metal
        principled.inputs['Roughness'].default_value = 0.35  # slight sheen

        if 'Specular IOR Level' in principled.inputs:
            principled.inputs['Specular IOR Level'].default_value = 0.5
        if 'IOR' in principled.inputs:
            principled.inputs['IOR'].default_value = 1.46  # plastic IOR

        # If there's subsurface scattering, enable for plastic translucency
        if 'Subsurface Weight' in principled.inputs:
            principled.inputs['Subsurface Weight'].default_value = 0.05
            principled.inputs['Subsurface Radius'].default_value = (0.8, 0.6, 0.4)
        elif 'Subsurface' in principled.inputs:
            principled.inputs['Subsurface'].default_value = 0.05

        print("  Base Color: bone white (0.85, 0.82, 0.78)")
        print("  Metallic: 0.0 (plastic)")
        print("  Roughness: 0.35")
    else:
        print("ERROR: No Principled BSDF found in material")

    # Find emission node and change cyan to warm white glow
    for n in nodes:
        if n.type == 'EMISSION':
            n.inputs['Color'].default_value = (1.0, 0.95, 0.85, 1.0)  # warm white glow
            print("  Emission: warm white (1.0, 0.95, 0.85)")

# --- 2. Set camera ortho_scale ---
cam = bpy.data.objects.get("Camera")
if cam and cam.data:
    cam.data.ortho_scale = 8.5
    print(f"  Camera ortho_scale: 8.5")

# --- 3. Set render resolution for UI half-scale ---
scene = bpy.context.scene
scene.render.resolution_x = 269
scene.render.resolution_y = 69
scene.render.image_settings.color_mode = 'RGBA'
scene.render.film_transparent = True
scene.cycles.samples = 64
scene.cycles.use_denoising = True
print(f"  Render: {scene.render.resolution_x}x{scene.render.resolution_y}")

# --- 4. Render 128 frames ---
out_dir = r"C:\Users\hooki\Trench\trench-plugin\assets\thumbwheel"
if not os.path.exists(out_dir):
    os.makedirs(out_dir)

scene.frame_start = 1
scene.frame_end = 128

print(f"\nRendering 128 frames to {out_dir}...")
for i in range(1, 129):
    scene.frame_set(i)
    scene.render.filepath = os.path.join(out_dir, f'frame_{i:04d}.png')
    bpy.ops.render.render(write_still=True)
    if i % 16 == 0:
        print(f"  Frame {i}/128")

print("Done.")
