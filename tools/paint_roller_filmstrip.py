"""
Paint TRENCH roller filmstrips at native resolution.

Heritage reference: E-mu BITMAP4330 — 128 frames, ~86×16 per frame.
Structure: cylindrical drum with knurl ridges above/below, emissive
circumferential band at the equator. Band wraps around cylinder surface —
foreshortens at edges, disappears behind. Two bands 180° apart so there's
always a position indicator visible.

Pixel-to-angle mapping uses arcsin for correct cylindrical projection.
"""

import math
import numpy as np
from PIL import Image

FRAME_W = 123
FRAME_H = 28
FRAMES = 128

# Vertical zone layout (rows)
SHOULDER_TOP = 2      # housing shadow
KNURL_TOP = 9         # upper knurl surface
BAND_ZONE = 8         # center emissive groove
KNURL_BOT = 7         # lower knurl surface
SHOULDER_BOT = 2      # housing shadow
# total = 2 + 9 + 8 + 7 + 2 = 28

# Cylinder
RIDGE_COUNT = 52       # ridges around full circumference
RIDGE_DEPTH = 0.30     # groove depth on knurl surface
BAND_HALF_ANGLE = 0.55 # band half-width in radians (~31°, so full band ~63°)

# Colors
KNURL_BASE = np.array([72, 66, 58], dtype=float)
KNURL_SPEC = np.array([165, 155, 140], dtype=float)
GROOVE_DARK = np.array([18, 16, 14], dtype=float)

# Primary band: TRENCH cyan
BAND1_CORE = np.array([160, 255, 255], dtype=float)
BAND1_MID  = np.array([45, 195, 200], dtype=float)
BAND1_EDGE = np.array([15, 80, 85], dtype=float)

# Secondary band (180° opposite): dimmer warm accent
BAND2_CORE = np.array([200, 140, 170], dtype=float)
BAND2_MID  = np.array([120, 55, 80], dtype=float)
BAND2_EDGE = np.array([50, 20, 35], dtype=float)

SHOULDER_COLOR = np.array([30, 28, 26], dtype=float)


def clamp8(v):
    return max(0, min(255, int(v + 0.5)))


def angle_dist(a, b):
    """Shortest angular distance on a circle."""
    d = (a - b) % (2 * math.pi)
    if d > math.pi:
        d = 2 * math.pi - d
    return d


def pixel_to_angle(x):
    """Map pixel x-position to angle on cylinder surface.
    Uses arcsin for correct cylindrical projection.
    Returns angle in [-pi/2, +pi/2]. Center pixel = 0."""
    xn = (x + 0.5) / FRAME_W  # 0..1
    xn = max(0.001, min(0.999, xn))  # clamp to avoid domain error
    return math.asin(2.0 * xn - 1.0)


def paint_frame(frame_idx):
    img = np.zeros((FRAME_H, FRAME_W, 4), dtype=np.uint8)

    # Drum rotation: full 360° over 128 frames
    rotation = (frame_idx / FRAMES) * 2.0 * math.pi

    # Band 1 is at angle 0 on the drum surface
    # Band 2 is at angle pi (opposite side)
    # After rotation, band 1 world angle = rotation, band 2 = rotation + pi

    # Specular highlight position (fixed in world, moves across drum as it rotates)
    spec_angle = 0.15  # slightly off-center for interest

    # Zone boundaries (row indices)
    y_knurl_top_start = SHOULDER_TOP
    y_knurl_top_end = SHOULDER_TOP + KNURL_TOP
    y_band_start = y_knurl_top_end
    y_band_end = y_band_start + BAND_ZONE
    y_knurl_bot_start = y_band_end
    y_knurl_bot_end = y_knurl_bot_start + KNURL_BOT

    for y in range(FRAME_H):
        # Vertical cylinder shading for cross-section
        # Map y within the drum body (excluding shoulders) to cylinder angle
        if y < SHOULDER_TOP:
            # Top shoulder: dark housing
            for x in range(FRAME_W):
                fade = (y + 1) / (SHOULDER_TOP + 1)
                c = SHOULDER_COLOR * fade
                img[y, x] = [clamp8(c[0]), clamp8(c[1]), clamp8(c[2]), int(200 * fade)]
            continue
        elif y >= FRAME_H - SHOULDER_BOT:
            # Bottom shoulder
            rows_from_bottom = FRAME_H - 1 - y
            fade = (rows_from_bottom + 1) / (SHOULDER_BOT + 1)
            for x in range(FRAME_W):
                c = SHOULDER_COLOR * fade
                img[y, x] = [clamp8(c[0]), clamp8(c[1]), clamp8(c[2]), int(200 * fade)]
            continue

        # Drum body rows: compute vertical shading
        drum_y = (y - SHOULDER_TOP) / (FRAME_H - SHOULDER_TOP - SHOULDER_BOT - 1)
        # Cylinder cross-section: bright at top (light from above), darker at bottom
        vert_shade = 0.45 + 0.55 * (1.0 - drum_y) ** 0.7
        # Specular band near top of drum
        vert_spec = max(0.0, 1.0 - abs(drum_y - 0.18) / 0.12) ** 1.8

        # Determine zone
        in_band_zone = y_band_start <= y < y_band_end
        # Transition rows at zone boundaries
        band_zone_t = 0.0
        if in_band_zone:
            band_zone_t = 1.0
            # Soft entry/exit at zone edges
            if y == y_band_start:
                band_zone_t = 0.6
            elif y == y_band_end - 1:
                band_zone_t = 0.5
            elif y == y_band_start + 1:
                band_zone_t = 0.85

        for x in range(FRAME_W):
            # Map pixel to angle on cylinder
            phi = pixel_to_angle(x)  # -pi/2 to +pi/2

            # Surface normal dot view direction (Lambertian)
            n_dot_v = math.cos(phi)

            # World angle of this surface point
            world_phi = phi + rotation

            # Ridge pattern (on the drum surface, wraps with rotation)
            ridge_phase = world_phi * RIDGE_COUNT
            ridge_cos = math.cos(ridge_phase)
            groove = max(0.0, -ridge_cos)
            peak = max(0.0, ridge_cos)

            # Specular highlight (broad, fixed in view space)
            spec_dist = abs(phi - spec_angle)
            spec_intensity = max(0.0, 1.0 - spec_dist / 0.6) ** 2.5

            if in_band_zone:
                # Check distance to each band
                d1 = angle_dist(world_phi, 0.0)       # band 1 at angle 0
                d2 = angle_dist(world_phi, math.pi)    # band 2 at angle pi

                # Band 1 (cyan)
                if d1 < BAND_HALF_ANGLE:
                    t = 1.0 - (d1 / BAND_HALF_ANGLE)
                    # Smooth hermite falloff
                    t = t * t * (3.0 - 2.0 * t)
                    if t > 0.7:
                        color = BAND1_MID + (BAND1_CORE - BAND1_MID) * ((t - 0.7) / 0.3)
                    else:
                        color = BAND1_EDGE + (BAND1_MID - BAND1_EDGE) * (t / 0.7)
                    # Emissive: partially ignores shading but not completely
                    emit = 0.4 + 0.6 * n_dot_v
                    color = color * emit * band_zone_t
                    # Ridges visible as subtle shadow on the band
                    color = color * (1.0 - 0.06 * groove)

                elif d2 < BAND_HALF_ANGLE:
                    t = 1.0 - (d2 / BAND_HALF_ANGLE)
                    t = t * t * (3.0 - 2.0 * t)
                    if t > 0.7:
                        color = BAND2_MID + (BAND2_CORE - BAND2_MID) * ((t - 0.7) / 0.3)
                    else:
                        color = BAND2_EDGE + (BAND2_MID - BAND2_EDGE) * (t / 0.7)
                    emit = 0.4 + 0.6 * n_dot_v
                    color = color * emit * band_zone_t
                    color = color * (1.0 - 0.06 * groove)

                else:
                    # Dark groove between bands
                    color = GROOVE_DARK * n_dot_v * vert_shade
                    # Mix with knurl at zone edges
                    if band_zone_t < 1.0:
                        knurl_c = KNURL_BASE * n_dot_v * vert_shade * (1.0 - RIDGE_DEPTH * groove)
                        knurl_c = knurl_c + KNURL_BASE * 0.3 * peak * n_dot_v
                        color = color * band_zone_t + knurl_c * (1.0 - band_zone_t)

                # Glow bleed from band into nearby rows
                if not in_band_zone:
                    pass  # handled below
            else:
                # Knurl zone
                ridge_factor = 1.0 - RIDGE_DEPTH * groove
                color = KNURL_BASE * n_dot_v * vert_shade * ridge_factor
                # Ridge peak highlight
                color = color + KNURL_BASE * 0.35 * peak * n_dot_v * vert_shade

                # Glow bleed from band (check proximity to band zone)
                glow_dist = 0
                if y < y_band_start:
                    glow_dist = y_band_start - y
                elif y >= y_band_end:
                    glow_dist = y - y_band_end + 1
                if glow_dist > 0 and glow_dist <= 3:
                    # Check if a band is facing us
                    d1 = angle_dist(world_phi, 0.0)
                    d2 = angle_dist(world_phi, math.pi)
                    glow_fade = (1.0 - glow_dist / 4.0) ** 2
                    if d1 < BAND_HALF_ANGLE * 0.8:
                        bt = 1.0 - d1 / (BAND_HALF_ANGLE * 0.8)
                        glow = BAND1_MID * bt * glow_fade * 0.25 * n_dot_v
                        color = color + glow
                    elif d2 < BAND_HALF_ANGLE * 0.8:
                        bt = 1.0 - d2 / (BAND_HALF_ANGLE * 0.8)
                        glow = BAND2_MID * bt * glow_fade * 0.15 * n_dot_v
                        color = color + glow

            # Specular on knurl (less on band zone)
            if not in_band_zone:
                color = color + KNURL_SPEC * spec_intensity * vert_spec * n_dot_v
            else:
                color = color + KNURL_SPEC * spec_intensity * vert_spec * n_dot_v * 0.15

            # Alpha
            a = 255

            img[y, x] = [clamp8(color[0]), clamp8(color[1]), clamp8(color[2]), a]

    return img


def main():
    import os
    import time

    out_dir = os.path.join(
        os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
        "trench-plugin", "assets",
    )
    thumb_dir = os.path.join(out_dir, "thumbwheel")
    os.makedirs(thumb_dir, exist_ok=True)

    print(f"Painting {FRAMES} frames at {FRAME_W}x{FRAME_H}...")
    t0 = time.time()

    frames = []
    for i in range(FRAMES):
        frame = paint_frame(i)
        frames.append(frame)
        if (i + 1) % 32 == 0:
            print(f"  frame {i+1}/{FRAMES}")

    # Stitch strip
    strip_w = FRAME_W * FRAMES
    strip = np.zeros((FRAME_H, strip_w, 4), dtype=np.uint8)
    for i, frame in enumerate(frames):
        strip[:, i * FRAME_W:(i + 1) * FRAME_W, :] = frame

    strip_path = os.path.join(out_dir, "thumbwheel_strip.png")
    Image.fromarray(strip).save(strip_path)
    print(f"Strip: {strip_path} ({strip_w}x{FRAME_H})")

    # Individual frames
    for i, frame in enumerate(frames):
        fpath = os.path.join(thumb_dir, f"frame_{i+1:04d}.png")
        Image.fromarray(frame).save(fpath)

    # 8x previews
    for i in [0, 15, 31, 47, 63, 79, 95, 111, 127]:
        frame_img = Image.fromarray(frames[i])
        up = frame_img.resize((FRAME_W * 8, FRAME_H * 8), Image.NEAREST)
        up_path = os.path.join(out_dir, f"thumbwheel_frame_{i+1:03d}_x8.png")
        up.save(up_path)
        print(f"Preview: {up_path}")

    elapsed = time.time() - t0
    print(f"Done in {elapsed:.1f}s")


if __name__ == "__main__":
    main()
