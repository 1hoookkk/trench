"""Stdlib-only PNG plotter test — is iOS rendering this?

Writes a minimal RGB PNG encoding a frequency-response curve. Uses
only `struct` and `zlib` from the Python standard library, so it
runs in any clean venv.
"""
from __future__ import annotations

import json
import math
import struct
import subprocess
import zlib
from pathlib import Path

# ---------------------------------------------------------------------------
# Minimal PNG encoder (RGB, 8-bit, filter=none)
# ---------------------------------------------------------------------------

def _chunk(tag: bytes, data: bytes) -> bytes:
    length = struct.pack(">I", len(data))
    body = tag + data
    crc = struct.pack(">I", zlib.crc32(body) & 0xFFFFFFFF)
    return length + body + crc


def write_png(path: Path, pixels: list, width: int, height: int) -> None:
    """pixels: flat bytearray of width*height*3 RGB bytes."""
    raw = bytearray()
    stride = width * 3
    for y in range(height):
        raw.append(0)  # filter type: none
        start = y * stride
        raw.extend(pixels[start : start + stride])
    compressed = zlib.compress(bytes(raw), 9)

    out = bytearray(b"\x89PNG\r\n\x1a\n")
    # IHDR: width, height, bit_depth=8, color_type=2 (RGB),
    #       compression=0, filter=0, interlace=0
    ihdr = struct.pack(">IIBBBBB", width, height, 8, 2, 0, 0, 0)
    out.extend(_chunk(b"IHDR", ihdr))
    out.extend(_chunk(b"IDAT", compressed))
    out.extend(_chunk(b"IEND", b""))
    path.write_bytes(bytes(out))


# ---------------------------------------------------------------------------
# Plotting primitives
# ---------------------------------------------------------------------------

class Canvas:
    def __init__(self, w: int, h: int, bg=(26, 26, 32)):
        self.w = w
        self.h = h
        self.px = bytearray(w * h * 3)
        for i in range(w * h):
            self.px[i * 3] = bg[0]
            self.px[i * 3 + 1] = bg[1]
            self.px[i * 3 + 2] = bg[2]

    def set(self, x: int, y: int, rgb: tuple) -> None:
        if 0 <= x < self.w and 0 <= y < self.h:
            i = (y * self.w + x) * 3
            self.px[i] = rgb[0]
            self.px[i + 1] = rgb[1]
            self.px[i + 2] = rgb[2]

    def hline(self, x0: int, x1: int, y: int, rgb: tuple) -> None:
        if x1 < x0:
            x0, x1 = x1, x0
        for x in range(x0, x1 + 1):
            self.set(x, y, rgb)

    def vline(self, x: int, y0: int, y1: int, rgb: tuple) -> None:
        if y1 < y0:
            y0, y1 = y1, y0
        for y in range(y0, y1 + 1):
            self.set(x, y, rgb)

    def line(self, x0: int, y0: int, x1: int, y1: int, rgb: tuple) -> None:
        # Bresenham
        dx = abs(x1 - x0)
        dy = -abs(y1 - y0)
        sx = 1 if x0 < x1 else -1
        sy = 1 if y0 < y1 else -1
        err = dx + dy
        while True:
            self.set(x0, y0, rgb)
            if x0 == x1 and y0 == y1:
                break
            e2 = 2 * err
            if e2 >= dy:
                err += dy
                x0 += sx
            if e2 <= dx:
                err += dx
                y0 += sy

    def thick_line(self, x0: int, y0: int, x1: int, y1: int, rgb: tuple) -> None:
        self.line(x0, y0, x1, y1, rgb)
        self.line(x0, y0 + 1, x1, y1 + 1, rgb)
        self.line(x0 + 1, y0, x1 + 1, y1, rgb)

    def filled_disc(self, cx: int, cy: int, r: int, rgb: tuple) -> None:
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                if dx * dx + dy * dy <= r * r:
                    self.set(cx + dx, cy + dy, rgb)


# ---------------------------------------------------------------------------
# DSP
# ---------------------------------------------------------------------------

SR = 39062.5

def biquad_mag(c, omega):
    cos1, sin1 = math.cos(omega), math.sin(omega)
    cos2, sin2 = math.cos(2 * omega), math.sin(2 * omega)
    num_re = c[0] + c[1] * cos1 + c[2] * cos2
    num_im = -c[1] * sin1 - c[2] * sin2
    den_re = 1.0 + c[3] * cos1 + c[4] * cos2
    den_im = -c[3] * sin1 - c[4] * sin2
    num = math.sqrt(num_re ** 2 + num_im ** 2)
    den = math.sqrt(den_re ** 2 + den_im ** 2)
    return num / den if den > 1e-18 else 1e6


def cascade_db(stages, boost, f):
    omega = 2 * math.pi * f / SR
    m = boost
    for s in stages:
        m *= biquad_mag(s, omega)
    return 20 * math.log10(m) if m > 1e-10 else -200.0


# ---------------------------------------------------------------------------
# Plot P2K M100/Q100 as a 640x360 PNG
# ---------------------------------------------------------------------------

def main():
    # Load P2K_013 from git history. The file was deleted in commit
    # f146fb3 ("drop pre-pill cartridge archives") so we extract via
    # `git show` rather than keeping a recovered copy in-tree.
    repo = Path(__file__).resolve().parents[1]
    result = subprocess.run(
        ["git", "-C", str(repo), "show", "f146fb3^:cartridges/p2k/P2k_013.json"],
        capture_output=True, text=True, check=True,
    )
    p2k = json.loads(result.stdout)

    # Pick M100_Q100
    kf = next(k for k in p2k["keyframes"] if k["label"] == "M100_Q100")
    stages = [[s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]] for s in kf["stages"]]
    boost = kf.get("boost", 1.0)

    # Canvas
    W, H = 640, 360
    cv = Canvas(W, H, bg=(26, 26, 32))

    # Plot area
    L, R, T, B = 60, 20, 40, 40
    pw = W - L - R
    ph = H - T - B
    db_min, db_max = -30, 45
    f_min, f_max = 20.0, 18000.0

    def x_for(f):
        return L + int((math.log10(f) - math.log10(f_min)) / (math.log10(f_max) - math.log10(f_min)) * pw)

    def y_for(db):
        clamped = max(db_min, min(db_max, db))
        return T + int((db_max - clamped) / (db_max - db_min) * ph)

    # Plot area background
    for y in range(T, T + ph + 1):
        cv.hline(L, L + pw, y, (14, 14, 20))

    # Grid
    grid = (42, 42, 53)
    for f in [50, 100, 200, 500, 1000, 2000, 5000, 10000]:
        cv.vline(x_for(f), T, T + ph, grid)
    for db in range(db_min, db_max + 1, 10):
        color = (85, 85, 85) if db == 0 else grid
        cv.hline(L, L + pw, y_for(db), color)

    # Plot area border
    border = (58, 58, 70)
    cv.hline(L, L + pw, T, border)
    cv.hline(L, L + pw, T + ph, border)
    cv.vline(L, T, T + ph, border)
    cv.vline(L + pw, T, T + ph, border)

    # Curve: sample 512 log-spaced frequencies
    curve_color = (255, 184, 74)  # amber
    N = 512
    prev = None
    for i in range(N):
        f = 10 ** (math.log10(f_min) + (math.log10(f_max) - math.log10(f_min)) * i / (N - 1))
        db = cascade_db(stages, boost, f)
        x = x_for(f)
        y = y_for(db)
        if prev is not None:
            cv.thick_line(prev[0], prev[1], x, y, curve_color)
        prev = (x, y)

    # Find and mark peaks
    peak_color = (255, 255, 255)
    peaks_hz = [194, 1509, 2139, 2411, 4399, 8988]
    for f in peaks_hz:
        db = cascade_db(stages, boost, f)
        cv.filled_disc(x_for(f), y_for(db), 3, peak_color)

    # Write out
    out_path = repo / "plots" / "p2k_013_m100_q100.png"
    out_path.parent.mkdir(exist_ok=True)
    write_png(out_path, cv.px, W, H)
    print(f"wrote {out_path}")
    print(f"size: {out_path.stat().st_size} bytes")


if __name__ == "__main__":
    main()
