#!/usr/bin/env python3
"""
UI Image Measurement Tool

Usage:
    python tools/measure_ui.py <image_path>...
    python tools/measure_ui.py *.png
    python tools/measure_ui.py docs/identity/*.png
"""

import sys
from pathlib import Path
from PIL import Image
import os


def measure_image(path: Path) -> dict:
    """Measure a single UI image."""
    try:
        with Image.open(path) as img:
            width, height = img.size
            mode = img.mode
            format_str = img.format or "UNKNOWN"
            
            # File size
            file_size = path.stat().st_size
            
            # Get basic stats if possible
            stats = {}
            try:
                if img.mode in ('RGB', 'RGBA', 'L'):
                    stat = img.getpixel((0, 0))
                    stats['has_transparency'] = img.mode == 'RGBA'
            except:
                pass
            
            return {
                'path': str(path),
                'width': width,
                'height': height,
                'mode': mode,
                'format': format_str,
                'file_size_bytes': file_size,
                'file_size_kb': file_size / 1024,
                'aspect_ratio': width / height if height > 0 else 0,
                'is_powers_of_two': ((width & (width - 1)) == 0 and (height & (height - 1)) == 0),
                'dpi': img.info.get('dpi', (72, 72)),
            }
    except Exception as e:
        return {
            'path': str(path),
            'error': str(e),
        }


def format_measurement(m: dict) -> str:
    """Format measurement for display."""
    if 'error' in m:
        return f"  ✗ {m['path']}: ERROR - {m['error']}"
    
    size_str = f"{m['width']}x{m['height']}"
    size_kb = f"{m['file_size_kb']:.1f}KB"
    ar = f"{m['aspect_ratio']:.2f}"
    pow2 = "✓" if m['is_powers_of_two'] else "✗"
    
    return f"  {m['path']}\n    size: {size_str} | {size_kb} | {m['mode']} | AR: {ar} | POW2: {pow2}"


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return
    
    # Collect all image paths
    paths = []
    for arg in sys.argv[1:]:
        p = Path(arg)
        if p.is_file():
            paths.append(p)
        elif p.is_dir():
            paths.extend(p.glob("*.png"))
            paths.extend(p.glob("*.jpg"))
            paths.extend(p.glob("*.jpeg"))
        else:
            # Try glob pattern
            parent = Path(".")
            paths.extend(parent.glob(arg))
    
    if not paths:
        print("No images found.")
        return
    
    print(f"Measuring {len(paths)} image(s)...\n")
    
    measurements = []
    for p in sorted(set(paths)):
        m = measure_image(p)
        measurements.append(m)
    
    for m in measurements:
        print(format_measurement(m))
    
    # Summary
    valid = [m for m in measurements if 'error' not in m]
    if len(valid) > 1:
        total_size = sum(m['file_size_bytes'] for m in valid)
        print(f"\n--- Summary ---")
        print(f"  {len(valid)} images, {total_size / 1024:.1f}KB total")


if __name__ == "__main__":
    main()