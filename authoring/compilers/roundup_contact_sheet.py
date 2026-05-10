import json
import math
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as signal

SR = 39062.5
IR_LEN = 2048
F_MIN = 40.0
F_MAX = 20000.0
DB_MIN = -60.0
DB_MAX = 60.0

BG = "#0b1020"
GRID = "#1e293b"
TICK = "#64748b"
TITLE = "#cbd5f5"

LAYOUT = (
    ("M0_Q0", "#7dd3fc"),
    ("M100_Q0", "#fb7185"),
)

def run_impulse(stages, length):
    x = np.zeros(length)
    x[0] = 1.0
    y = x.copy()
    for s in stages:
        # Standard filter form in TRENCH:
        # y[n] = c0*x[n] + s1[n-1]
        # s1[n] = c1*x[n] - c3*y[n] + s2[n-1]
        # s2[n] = c2*x[n] - c4*y[n]
        # This is matches tools/plot_magnitude.py math:
        # H(z) = (c0 + c1 z^-1 + c2 z^-2) / (1 + c3 z^-1 + c4 z^-2)
        b = [s[0], s[1], s[2]]
        a = [1.0, s[3], s[4]]
        y = signal.lfilter(b, a, y)
    return y

def get_curve(kf):
    stages = []
    # Support both compiled (c0..c4) and raw (a1..val)
    if "c0" in kf["stages"][0]:
        stages = [[s["c0"], s["c1"], s["c2"], s["c3"], s["c4"]] for s in kf["stages"]]
    else:
        # Raw -> Encoded conversion (simplified)
        for s in kf["stages"]:
            a1, r, v1, v2, v3 = s["a1"], s["r"], s["val1"], s["val2"], s["val3"]
            a2 = r * r
            b0 = 1.0 + v1
            b1 = a1 + v2
            b2 = a2 - v3
            stages.append([b0, b1, b2, a1, a2])
            
    ir = run_impulse(stages, IR_LEN)
    boost = float(kf.get("boost", 1.0))
    w, h = signal.freqz(ir, worN=1024, fs=SR)
    mag = np.abs(h) * boost
    mag = np.nan_to_num(mag, nan=1e-12, posinf=1e6, neginf=1e-12)
    mag_db = 20 * np.log10(np.maximum(mag, 1e-12))
    return w, mag_db

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("dir", type=Path)
    parser.add_argument("--out-prefix", default="roundup_sheet")
    args = parser.parse_args()
    
    files = sorted([f for f in args.dir.glob("*.json") if f.is_file()])
    print(f"Found {len(files)} files")
    
    batch_size = 12
    num_batches = math.ceil(len(files) / batch_size)
    
    for b in range(num_batches):
        batch_files = files[b*batch_size : (b+1)*batch_size]
        cols = 3
        rows = math.ceil(len(batch_files) / cols)
        
        fig, axes = plt.subplots(rows, cols, figsize=(cols*5, rows*4), dpi=120, facecolor=BG)
        axes_arr = np.atleast_1d(axes).flatten()
        
        for ax, f_path in zip(axes_arr, batch_files):
            with open(f_path) as f:
                doc = json.load(f)
            
            kfs = {kf["label"]: kf for kf in doc.get("keyframes", [])}
            if not kfs and "corners" in doc:
                kfs = {k: {"stages": v["stages"], "label": k} for k, v in doc["corners"].items()}
            
            for label, color in LAYOUT:
                if label in kfs:
                    w, mag_db = get_curve(kfs[label])
                    ax.semilogx(w, mag_db, color=color, linewidth=2.0, alpha=0.9)
            
            ax.set_title(f_path.stem, color=TITLE, fontsize=10, pad=2)
            ax.set_facecolor(BG)
            ax.set_ylim(DB_MIN, DB_MAX)
            ax.set_xlim(F_MIN, F_MAX)
            ax.grid(True, which="both", color=GRID, alpha=0.4, linestyle=":")
            ax.tick_params(colors=TICK, labelsize=7, length=2)
            for sp in ax.spines.values():
                sp.set_color(GRID)
            # Remove axis labels for contact sheet density
            ax.set_xlabel("")
            ax.set_ylabel("")
        
        for ax in axes_arr[len(batch_files):]:
            ax.axis("off")
            
        fig.tight_layout()
        out_path = Path(f"{args.out_prefix}_{b:02d}.png")
        fig.savefig(out_path, facecolor=BG)
        plt.close(fig)
        print(f"Wrote {out_path}")

if __name__ == "__main__":
    main()
