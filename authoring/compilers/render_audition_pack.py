import subprocess
import os
import sys
import numpy as np
import wave
from pathlib import Path

def generate_saw(path, duration=2.0, freq=100.0, sr=44100):
    t = np.linspace(0, duration, int(sr * duration), False)
    saw = 2.0 * (freq * t - np.floor(0.5 + freq * t))
    
    with wave.open(str(path), 'wb') as wav:
        wav.setnchannels(1)
        wav.setsampwidth(4) # 32-bit
        wav.setframerate(sr)
        # Manually pack as f32 if wave module doesn't support it directly
        # But wait, wave.py usually handles PCM. For IEEE float, we need a different approach or a library.
        # I'll use a simple manual pack.
        import struct
        data = b"".join(struct.pack('f', s) for s in saw)
        # Note: wave module might write the wrong header if we just give it bytes.
        # Actually, let's just use the renderer's expected format.
        # The renderer expects RIFF/WAVE with fmt chunk and data chunk.
    
    # I'll rewrite the generators to use a more robust way to write f32 wavy
    write_f32_wav(path, saw, sr)

def write_f32_wav(path, data, sr):
    import struct
    content = bytearray()
    content.extend(b"RIFF")
    content.extend(struct.pack("<I", 36 + len(data)*4))
    content.extend(b"WAVE")
    content.extend(b"fmt ")
    content.extend(struct.pack("<I", 16))
    content.extend(struct.pack("<H", 3)) # IEEE float
    content.extend(struct.pack("<H", 1)) # Mono
    content.extend(struct.pack("<I", sr))
    content.extend(struct.pack("<I", sr * 4))
    content.extend(struct.pack("<H", 4))
    content.extend(struct.pack("<H", 32))
    content.extend(b"data")
    content.extend(struct.pack("<I", len(data)*4))
    for s in data:
        content.extend(struct.pack("<f", float(s)))
    with open(path, "wb") as f:
        f.write(content)

def generate_noise(path, duration=2.0, sr=44100):
    noise = np.random.uniform(-1, 1, int(sr * duration))
    write_f32_wav(path, noise, sr)

def generate_sine(path, duration=2.0, freq=50.0, sr=44100):
    t = np.linspace(0, duration, int(sr * duration), False)
    sine = np.sin(2 * np.pi * freq * t)
    write_f32_wav(path, sine, sr)

def render_pack():
    base_dir = Path(r"c:\Users\hooki\Trench")
    roundup_dir = base_dir / "roundup"
    out_dir = base_dir / "auditions" / "four_doors"
    out_dir.mkdir(parents=True, exist_ok=True)
    
    # 1. Setup Sources
    src_dir = out_dir / "sources"
    src_dir.mkdir(exist_ok=True)
    
    saw_path = src_dir / "saw.wav"
    noise_path = src_dir / "noise.wav"
    sine50_path = src_dir / "sine50.wav"
    vocal_src = Path(r"c:\Users\hooki\trenchwork\datasets\wav_corpus\vocals_vocalset\vocal_sing_f1_long_straight_a.wav")
    vocal_dest = src_dir / "vocal.wav"
    
    if not saw_path.exists(): generate_saw(saw_path)
    if not noise_path.exists(): generate_noise(noise_path)
    if not sine50_path.exists(): generate_sine(sine50_path)
    # Copy vocal if possible, or just use path
    
    # 2. Define Jobs
    # (Body, Source, Note)
    jobs = [
        # Small Talk
        ("Vowel_Ah_to_Eh", saw_path, "Small_Talk_AhEe"),
        ("Vowel_Ah_to_Eh", str(vocal_src), "Small_Talk_AhEe"),
        
        # Speaker Knockerz
        ("Contra_octave_chain", sine50_path, "Speaker_Knockerz"),
        ("Contra_octave_chain", saw_path, "Speaker_Knockerz"),
        
        # Aluminum Siding
        ("Iron_Float", noise_path, "Aluminum_Siding"),
        ("Iron_Float", saw_path, "Aluminum_Siding"),
        
        # Cul-De-Sac
        ("Descending_Ladder", saw_path, "Cul_De_Sac"),
        ("Descending_Ladder", noise_path, "Cul_De_Sac")
    ]
    
    render_exe = base_dir / "trench-juce" / "plugin" / "build" / "TRENCH_OfflineRender_artefacts" / "Release" / "TRENCH_OfflineRender.exe"
    if not render_exe.exists():
        # Try Debug if Release missing
        render_exe = base_dir / "trench-juce" / "plugin" / "build" / "TRENCH_OfflineRender_artefacts" / "Debug" / "TRENCH_OfflineRender.exe"
    
    print(f"Using renderer: {render_exe}")
    
    for body_name, src_path, label in jobs:
        body_path = roundup_dir / f"{body_name}.json"
        if not body_path.exists():
            print(f"DANGER: Body {body_name} not found in roundup!")
            continue
            
        src_name = Path(src_path).stem
        out_wav = out_dir / f"{label}_{src_name}.wav"
        
        print(f"Rendering {label} with {src_name}...")
        
        # Command: renderer.exe <cartridge_json> <input_wav> <output_wav> [morph_start] [morph_end] [q_val]
        # We'll do a 0->100 morph over 2 seconds
        # Command: renderer.exe <cartridge_json> <input.wav> <output.wav>
        cmd = [
            str(render_exe),
            str(body_path),
            str(src_path),
            str(out_wav)
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
        except Exception as e:
            print(f"Failed to render {label}: {e}")
            
    # 3. Generate Visuals
    plot_script = base_dir / "tools" / "plot_magnitude.py"
    for body_name, _, label in jobs:
        body_path = roundup_dir / f"{body_name}.json"
        # Plot Cascade
        out_plot = out_dir / f"{label}_magnitude.png"
        if not out_plot.exists():
            cmd = [sys.executable, str(plot_script), str(body_path), "--out", str(out_plot)]
            subprocess.run(cmd, capture_output=True)
            
    print(f"\nAudition pack ready in {out_dir}")

if __name__ == "__main__":
    render_pack()
