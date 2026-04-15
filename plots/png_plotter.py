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

def render_corner(repo: Path, p2k: dict, corner_label: str, out_name: str) -> None:
    """Render one P2K keyframe at the given corner as a 400x300 PNG
    in E-mu display style."""
    kf = next(k for k in p2k["keyframes"] if k["label"] == corner_label)
    stages = [[s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]] for s in kf["stages"]]
    boost = kf.get("boost", 1.0)
    _render_png(repo, stages, boost, out_name)


def _render_png(repo: Path, stages: list, boost: float, out_name: str) -> None:

    # Canvas — roughly matches the E-mu UI widget aspect ratio
    # (the plugin's filter display is ~250x180 pixels). Squarer layout
    # puts the 6 formant peaks visually where they sit in the UI.
    W, H = 400, 300
    # E-mu display color scheme: dark teal background, cyan grid,
    # bright cyan curve. Matches the Proteus 2000 plugin filter widget.
    BG        = (10, 28, 32)     # dark teal
    PLOT_BG   = (4, 18, 22)      # slightly darker inside the plot area
    GRID_MAJ  = (40, 90, 100)    # major gridlines (octave / 12 dB)
    GRID_MIN  = (20, 50, 58)     # minor gridlines (between octaves)
    BORDER    = (60, 130, 140)   # plot area border
    CURVE     = (106, 252, 216)  # bright cyan curve
    ZERO_AX   = (80, 180, 190)   # 0 dB axis highlight
    cv = Canvas(W, H, bg=BG)

    # Plot area
    L, R, T, B = 40, 15, 25, 25
    pw = W - L - R
    ph = H - T - B
    # Narrower dB range to match E-mu's plugin display. Real UI
    # clamps anything outside roughly +/- 24 dB so deep notches and
    # tall peaks both register at the edges without making the rest
    # of the curve look squashed. This is why the UI at M0/Q0
    # reads as a smooth lowpass even though the underlying math has
    # notches dipping to -60 dB — they just clip against the floor.
    db_min, db_max = -24, 24
    f_min, f_max = 50.0, 12000.0  # UI widget roughly spans this range

    def x_for(f):
        return L + int((math.log10(f) - math.log10(f_min)) / (math.log10(f_max) - math.log10(f_min)) * pw)

    def y_for(db):
        clamped = max(db_min, min(db_max, db))
        return T + int((db_max - clamped) / (db_max - db_min) * ph)

    # Plot area background
    for y in range(T, T + ph + 1):
        cv.hline(L, L + pw, y, PLOT_BG)

    # Logarithmic frequency grid — octave-spaced majors (every 2x),
    # with 1/3-octave minors between them (E-mu-style density).
    # Majors: powers of 2 from 50 Hz up to Nyquist.
    major_freqs = []
    f = 50.0
    while f <= f_max:
        major_freqs.append(f)
        f *= 2.0
    # Minors: 1.26x (1/3 octave) between each major
    minor_freqs = []
    for mf in major_freqs[:-1]:
        minor_freqs.append(mf * 1.2599)  # 2^(1/3)
        minor_freqs.append(mf * 1.5874)  # 2^(2/3)

    for mf in minor_freqs:
        if f_min <= mf <= f_max:
            cv.vline(x_for(mf), T, T + ph, GRID_MIN)
    for mf in major_freqs:
        if f_min <= mf <= f_max:
            cv.vline(x_for(mf), T, T + ph, GRID_MAJ)

    # dB grid — majors every 12 dB, minors every 6 dB, 0 dB
    # highlighted in the axis color so you can see the "rest"
    # level at a glance.
    for dbv in range(-24, 25, 6):
        if dbv == 0:
            color = ZERO_AX
        elif dbv % 12 == 0:
            color = GRID_MAJ
        else:
            color = GRID_MIN
        cv.hline(L, L + pw, y_for(dbv), color)

    # Plot area border
    cv.hline(L, L + pw, T, BORDER)
    cv.hline(L, L + pw, T + ph, BORDER)
    cv.vline(L, T, T + ph, BORDER)
    cv.vline(L + pw, T, T + ph, BORDER)

    # Curve: sample 4096 log-spaced frequencies. At 512 points we
    # were undersampling narrow high-Q resonators — adjacent samples
    # landed on opposite sides of a sharp peak and thick_line drew
    # jagged V-shapes where the true response is smooth. 4096 gives
    # ~1200 points per decade, enough that even 20 Hz bandwidth
    # peaks get 5+ samples across their main lobe.
    curve_color = CURVE  # bright cyan
    N = 4096
    prev = None
    for i in range(N):
        f = 10 ** (math.log10(f_min) + (math.log10(f_max) - math.log10(f_min)) * i / (N - 1))
        db = cascade_db(stages, boost, f)
        x = x_for(f)
        y = y_for(db)
        if prev is not None:
            cv.thick_line(prev[0], prev[1], x, y, curve_color)
        prev = (x, y)

    # No peak markers — E-mu's display shows the curve alone, so
    # leave it clean for direct visual comparison. If you need the
    # formant frequencies, the earlier console output has them.

    # Write out
    out_path = repo / "plots" / out_name
    out_path.parent.mkdir(exist_ok=True)
    write_png(out_path, cv.px, W, H)
    print(f"wrote {out_path.relative_to(repo)}  ({out_path.stat().st_size} bytes)")


def main():
    # All 33 P2K filters (P2k_000–P2k_032) live in git history under
    # cartridges/p2k/ — deleted in f146fb3 ("drop pre-pill cartridge
    # archives"). Extract each via `git show` and render all 4 corners.
    repo = Path(__file__).resolve().parents[1]

    corners = ["M0_Q0", "M100_Q0", "M0_Q100", "M100_Q100"]

    for n in range(33):
        num = f"{n:03d}"
        result = subprocess.run(
            ["git", "-C", str(repo), "show",
             f"f146fb3^:cartridges/p2k/P2k_{num}.json"],
            capture_output=True, text=True, check=True,
        )
        p2k = json.loads(result.stdout)
        for corner in corners:
            out_name = f"p2k_{num}_{corner.lower()}.png"
            render_corner(repo, p2k, corner, out_name)


if __name__ == "__main__":
    main()
