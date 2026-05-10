import os
import argparse
import json
import math
from PIL import Image, ImageDraw, ImageFilter, ImageOps

def lerp(a, b, t):
    return a + (b - a) * t

def create_rounded_rect_mask(size, radius):
    mask = Image.new('L', size, 0)
    draw = ImageDraw.Draw(mask)
    draw.rounded_rectangle((0, 0, size[0], size[1]), radius, fill=255)
    return mask

def build_rotor_strip():
    parser = argparse.ArgumentParser(description="TRENCH Rotor Sprite Strip Generator")
    parser.add_argument("--base", required=True, help="Path to inactive rotor PNG")
    parser.add_argument("--light-mask", help="Path to light mask PNG")
    parser.add_argument("--slot-mask", help="Path to slot mask PNG")
    parser.add_argument("--out-dir", required=True, help="Output directory for frames")
    parser.add_argument("--frames", type=int, default=128, help="Total frames")
    parser.add_argument("--strip", help="Output path for sprite strip")
    parser.add_argument("--contact", help="Output path for contact sheet")
    
    parser.add_argument("--width", type=int, default=300)
    parser.add_argument("--height", type=int, default=76)
    
    parser.add_argument("--cyan-strength", type=float, default=1.0)
    parser.add_argument("--glow-radius", type=int, default=12)
    parser.add_argument("--core-width", type=int, default=54)
    parser.add_argument("--core-height", type=int, default=12)
    parser.add_argument("--top-highlight-strength", type=float, default=0.35)
    parser.add_argument("--top-highlight-offset", type=int, default=8)
    parser.add_argument("--left-margin", type=int, default=28)
    parser.add_argument("--right-margin", type=int, default=272)
    
    args = parser.parse_args()

    if not os.path.exists(args.out_dir):
        os.makedirs(args.out_dir)

    # 1. Load Assets
    base_img = Image.open(args.base).convert("RGBA")
    if args.width and args.height:
        if base_img.size != (args.width, args.height):
            base_img = base_img.resize((args.width, args.height), Image.Resampling.LANCZOS)
    
    w, h = base_img.size
    
    light_mask = None
    if args.light_mask and os.path.exists(args.light_mask):
        light_mask = Image.open(args.light_mask).convert("L")
        if light_mask.size != (w, h):
            light_mask = light_mask.resize((w, h), Image.Resampling.LANCZOS)
    
    slot_mask = None
    if args.slot_mask and os.path.exists(args.slot_mask):
        slot_mask = Image.open(args.slot_mask).convert("L")
        if slot_mask.size != (w, h):
            slot_mask = slot_mask.resize((w, h), Image.Resampling.LANCZOS)

    # Fallback light mask logic
    if not light_mask:
        # Generate approximate mask: focus on mid-height, avoid brightest tooth tips
        light_mask = Image.new("L", (w, h), 0)
        draw = ImageDraw.Draw(light_mask)
        # Mid band area
        draw.rectangle((0, h//4, w, 3*h//4), fill=180)
        # Blur it
        light_mask = light_mask.filter(ImageFilter.GaussianBlur(radius=4))

    frames_images = []

    for f in range(args.frames):
        frame_canvas = base_img.copy()
        
        if 1 <= f <= 126:
            # Calculate position
            t = (f - 1) / 125.0
            x_center = lerp(args.left_margin, args.right_margin, t)
            
            # --- Cyan Glow Layer ---
            glow_layer = Image.new("RGBA", (w, h), (0, 0, 0, 0))
            glow_draw = ImageDraw.Draw(glow_layer)
            
            # Draw Core
            cw, ch = args.core_width, args.core_height
            cx1, cy1 = x_center - cw//2, (h//2) - ch//2
            cx2, cy2 = x_center + cw//2, (h//2) + ch//2
            
            # Core color #00E5E5
            core_color = (0, 229, 229, int(255 * args.cyan_strength))
            glow_draw.rounded_rectangle((cx1, cy1, cx2, cy2), radius=ch//2, fill=core_color)
            
            # Apply Blur for Glow
            glow_layer = glow_layer.filter(ImageFilter.GaussianBlur(radius=args.glow_radius))
            
            # Apply Light Mask (Occlusion)
            if light_mask:
                # Multiply glow alpha by light_mask
                glow_r, glow_g, glow_b, glow_a = glow_layer.split()
                new_a = ImageOps.procrustes = ImageChops_multiply_L(glow_a, light_mask) if 'ImageChops' in globals() else None
                # pillow ImageChops approach
                from PIL import ImageChops
                glow_a = ImageChops.multiply(glow_a, light_mask)
                glow_layer = Image.merge("RGBA", (glow_r, glow_g, glow_b, glow_a))

            # Composite Cyan Glow
            frame_canvas = Image.alpha_composite(frame_canvas, glow_layer)
            
            # --- Top Highlight Layer ---
            high_layer = Image.new("RGBA", (w, h), (0, 0, 0, 0))
            high_draw = ImageDraw.Draw(high_layer)
            
            hx = x_center + args.top_highlight_offset
            hw, hh = args.core_width * 0.8, 3
            hy = (h//3) # Upper third
            hx1, hy1 = hx - hw//2, hy - hh//2
            hx2, hy2 = hx + hw//2, hy + hh//2
            
            high_color = (230, 240, 255, int(255 * args.top_highlight_strength))
            high_draw.ellipse((hx1, hy1, hx2, hy2), fill=high_color)
            high_layer = high_layer.filter(ImageFilter.GaussianBlur(radius=4))
            
            frame_canvas = Image.alpha_composite(frame_canvas, high_layer)

        # Apply Slot Mask Clipping
        if slot_mask:
            r, g, b, a = frame_canvas.split()
            a = ImageChops.multiply(a, slot_mask)
            frame_canvas = Image.merge("RGBA", (r, g, b, a))

        # Save Frame
        frame_name = f"trench_rotor_{f:03d}.png"
        frame_path = os.path.join(args.out_dir, frame_name)
        frame_canvas.save(frame_path)
        frames_images.append(frame_canvas)

    # 2. Build Sprite Strip
    if args.strip:
        strip_img = Image.new("RGBA", (w * args.frames, h))
        for i, img in enumerate(frames_images):
            strip_img.paste(img, (i * w, 0))
        strip_img.save(args.strip)

    # 3. Build Contact Sheet
    if args.contact:
        cols = 16
        rows = math.ceil(args.frames / cols)
        contact_w, contact_h = w * cols, h * rows
        contact_img = Image.new("RGBA", (contact_w, contact_h), (0,0,0,0))
        for i, img in enumerate(frames_images):
            c = i % cols
            r = i // cols
            contact_img.paste(img, (c * w, r * h))
        contact_img.save(args.contact)

    # 4. Manifest
    manifest = {
        "frames": args.frames,
        "frame_width": w,
        "frame_height": h,
        "strip": os.path.basename(args.strip) if args.strip else None,
        "frame_0": "inactive",
        "frame_127": "inactive",
        "active_frames": [1, 126]
    }
    manifest_path = os.path.join(args.out_dir, "manifest.json")
    with open(manifest_path, "w") as f_man:
        json.dump(manifest, f_man, indent=2)

    print(f"Generated {args.frames} frames in {args.out_dir}")

if __name__ == "__main__":
    build_rotor_strip()
