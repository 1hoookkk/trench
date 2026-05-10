import bpy
import math

def build_glow_band():
    wheel = bpy.data.objects.get('Thumbwheel')
    if not wheel:
        print("Error: Could not find 'Thumbwheel' object.")
        return

    # 1. Spawn inner cylinder just underneath ribs
    # Assuming wheel radius is roughly 4.015m and Z is 1.5m
    bpy.ops.mesh.primitive_cylinder_add(
        vertices=128, 
        radius=4.01 - 0.04,  # Just inside outer geometry
        depth=1.45,          # Slightly narrower 
        location=wheel.location
    )
    inner_cyl = bpy.context.active_object
    inner_cyl.name = "InnerGlow"
    inner_cyl.rotation_euler = wheel.rotation_euler

    # 2. Parent it to thumbwheel so it native-spins
    inner_cyl.parent = wheel
    inner_cyl.matrix_parent_inverse = wheel.matrix_world.inverted()

    # 3. Create the Shader
    mat = bpy.data.materials.new(name="CyanBandShader")
    mat.use_nodes = True
    inner_cyl.data.materials.append(mat)

    nodes = mat.node_tree.nodes
    links = mat.node_tree.links
    for n in nodes:
        nodes.remove(n)

    mat_out = nodes.new('ShaderNodeOutputMaterial')
    mat_out.location = (1200, 0)
    
    emission = nodes.new('ShaderNodeEmission')
    emission.location = (900, 0)
    emission.inputs['Strength'].default_value = 8.0  # Hot emission inside the drum
    links.new(emission.outputs['Emission'], mat_out.inputs['Surface'])

    # UV Coordinates
    tex_coord = nodes.new('ShaderNodeTexCoord')
    tex_coord.location = (-1000, 0)
    
    sep_xyz = nodes.new('ShaderNodeSeparateXYZ')
    sep_xyz.location = (-800, 0)
    links.new(tex_coord.outputs['UV'], sep_xyz.inputs['Vector'])

    # --- U-AXIS LOGIC ---
    u_node = sep_xyz.outputs['X']
    
    # u - 0.58
    sub_u = nodes.new('ShaderNodeMath')
    sub_u.operation = 'SUBTRACT'
    sub_u.inputs[1].default_value = 0.58
    links.new(u_node, sub_u.inputs[0])
    
    # abs(u - 0.58)
    abs_u = nodes.new('ShaderNodeMath')
    abs_u.operation = 'ABSOLUTE'
    links.new(sub_u.outputs['Value'], abs_u.inputs[0])
    
    # 1 - abs(...)
    sub_1 = nodes.new('ShaderNodeMath')
    sub_1.operation = 'SUBTRACT'
    sub_1.inputs[0].default_value = 1.0
    links.new(abs_u.outputs['Value'], sub_1.inputs[1])
    
    # min(...)
    min_dist = nodes.new('ShaderNodeMath')
    min_dist.operation = 'MINIMUM'
    links.new(abs_u.outputs['Value'], min_dist.inputs[0])
    links.new(sub_1.outputs['Value'], min_dist.inputs[1])
    
    # Smoothstep band computation
    # width L = 0.24, s = 0.04
    # From Min = 0.12, From Max = 0.16. To Min = 1, To Max = 0 creates 1 - smoothstep out of the box
    band_range = nodes.new('ShaderNodeMapRange')
    band_range.interpolation_type = 'SMOOTHSTEP'
    band_range.inputs['From Min'].default_value = 0.12
    band_range.inputs['From Max'].default_value = 0.16
    band_range.inputs['To Min'].default_value = 1.0
    band_range.inputs['To Max'].default_value = 0.0
    links.new(min_dist.outputs['Value'], band_range.inputs['Value'])

    # --- V-AXIS LOGIC (VERTICAL MASK) ---
    v_node = sep_xyz.outputs['Y']
    v_ramp = nodes.new('ShaderNodeValToRGB')
    v_ramp.color_ramp.interpolation = 'B_SPLINE'
    v_ramp.color_ramp.elements[0].position = 0.20
    v_ramp.color_ramp.elements[0].color = (0, 0, 0, 1)
    
    el2 = v_ramp.color_ramp.elements.new(0.40)
    el2.color = (1, 1, 1, 1)
    
    el3 = v_ramp.color_ramp.elements.new(0.60)
    el3.color = (1, 1, 1, 1)
    
    v_ramp.color_ramp.elements[3].position = 0.80
    v_ramp.color_ramp.elements[3].color = (0, 0, 0, 1)
    
    links.new(v_node, v_ramp.inputs['Fac'])

    # Multiply Band x Mask
    mask_mult = nodes.new('ShaderNodeMath')
    mask_mult.operation = 'MULTIPLY'
    links.new(band_range.outputs['Result'], mask_mult.inputs[0])
    links.new(v_ramp.outputs['Color'], mask_mult.inputs[1])

    # --- COLOR MAPPING ---
    color_ramp = nodes.new('ShaderNodeValToRGB')
    color_ramp.color_ramp.interpolation = 'EASE'
    color_ramp.color_ramp.elements[0].position = 0.0
    color_ramp.color_ramp.elements[0].color = (0, 0, 0, 1)
    
    el4 = color_ramp.color_ramp.elements.new(0.3)
    el4.color = (0.01, 0.4, 0.5, 1) # Outer blue-green glow
    
    color_ramp.color_ramp.elements[2].position = 0.9
    color_ramp.color_ramp.elements[2].color = (0.2, 1.0, 1.0, 1) # Hot cyan core
    
    links.new(mask_mult.outputs['Value'], color_ramp.inputs['Fac'])
    links.new(color_ramp.outputs['Color'], emission.inputs['Color'])

    print("Success: Generated 'InnerGlow' emissive cylinder tied to drum phase!")

build_glow_band()
