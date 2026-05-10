#!/usr/bin/env python3
import asyncio
import json
import os
import tempfile
import sys
from pathlib import Path

# Add this directory to path to import compile_raw
tools_dir = Path(__file__).parent.absolute()
sys.path.append(str(tools_dir))

try:
    import websockets
except ImportError:
    print("[Forge] Error: 'websockets' library not found.")
    print("Please run: pip install websockets")
    sys.exit(1)

from compile_raw import compile_raw, CompileError

PROJECT_ROOT = tools_dir.parent
CARTRIDGE_PATH = PROJECT_ROOT / "active_forge_body.cartridge"

async def handle_forge(websocket):
    async for message in websocket:
        try:
            # 1. Parse JSON payload from UI
            payload = json.loads(message)
            name = payload.get("name", "unnamed_body")
            print(f"[Forge] Received update -> {name}")
            
            # 2. Compile to compiled-v1 (uint16 'w' arrays)
            try:
                compiled = compile_raw(payload)
            except CompileError as ce:
                print(f"[Forge] Compile Failed: {ce}")
                continue

            # 3. Atomic Write: temp file -> replace
            # This ensures the JUCE ArcSwap thread never reads a partial file.
            fd, temp_path = tempfile.mkstemp(dir=PROJECT_ROOT, suffix=".tmp")
            try:
                with os.fdopen(fd, 'w') as f:
                    json.dump(compiled, f, indent=2)
                
                os.replace(temp_path, CARTRIDGE_PATH)
                print(f"[Forge] Compiled -> Swapped {CARTRIDGE_PATH.name}")
            except Exception as e:
                if os.path.exists(temp_path):
                    os.remove(temp_path)
                raise e

        except json.JSONDecodeError:
            print("[Forge] Error: Received malformed JSON")
        except Exception as e:
            print(f"[Forge] Internal Error: {e}")

async def main():
    print(f"[Forge] Target: {CARTRIDGE_PATH}")
    async with websockets.serve(handle_forge, "localhost", 8765):
        print("[Forge] Server listening on ws://localhost:8765")
        await asyncio.Future()  # Run forever

if __name__ == "__main__":
    try:
        asyncio.run(main())
    except KeyboardInterrupt:
        print("\n[Forge] Server shutting down.")
