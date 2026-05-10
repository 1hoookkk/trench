#!/usr/bin/env python3
"""
TRENCH Body Baker - Test WebSocket Server
Simple echo server for testing the HTML body baker
"""
import asyncio
import json
import struct

async def handle_client(reader, writer):
    addr = writer.get_extra_info('peername')
    print(f"[TEST] Client connected: {addr}")
    
    buffer = b''
    
    try:
        while True:
            data = await reader.read(65536)
            if not data:
                break
            
            buffer += data
            
            # Simple WebSocket frame parser
            while len(buffer) >= 2:
                first = buffer[0]
                second = buffer[1]
                opcode = first & 0x0F
                masked = (second & 0x80) != 0
                payload_len = second & 0x7F
                header_len = 2 + (4 if masked else 0)
                
                if payload_len == 126 and len(buffer) >= 4:
                    payload_len = struct.unpack('>H', buffer[2:4])[0]
                    header_len = 4 + (4 if masked else 0)
                elif payload_len == 127 and len(buffer) >= 10:
                    payload_len = struct.unpack('>Q', buffer[2:10])[0]
                    header_len = 10 + (4 if masked else 0)
                
                if len(buffer) < header_len + payload_len:
                    break
                
                payload = buffer[header_len:header_len + payload_len]
                buffer = buffer[header_len + payload_len:]
                
                if opcode == 0x1:  # Text frame
                    try:
                        msg = json.loads(payload.decode('utf-8'))
                        print(f"[TEST] Body received: {msg.get('name', 'unnamed')}")
                        
                        # Count stages and corners
                        corners = msg.get('corners', {})
                        total_stages = sum(len(v.get('stages', [])) for v in corners.values())
                        print(f"[TEST] Corners: {len(corners)}, Stages: {total_stages}")
                        
                        # Echo response
                        response = json.dumps({
                            'status': 'received',
                            'name': msg.get('name'),
                            'stages_per_corner': len(list(corners.values())[0].get('stages', [])) if corners else 0
                        }).encode()
                        
                        # Send WebSocket frame
                        frame = bytearray()
                        frame.append(0x81)
                        if len(response) < 126:
                            frame.append(len(response))
                        else:
                            frame.extend([0x7E, *struct.pack('>H', len(response))])
                        
                        writer.write(bytes(frame) + response)
                        await writer.drain()
                        
                        print("[TEST] Response sent")
                        
                    except json.JSONDecodeError as e:
                        print(f"[TEST] Parse error: {e}")
                
                elif opcode == 0x8:
                    print(f"[TEST] Client closing: {addr}")
                    break
                    
    except Exception as e:
        print(f"[TEST] Error: {e}")
    finally:
        writer.close()
        await writer.wait_closed()
        print(f"[TEST] Client disconnected: {addr}")

async def main():
    print("=" * 50)
    print("TRENCH BODY BAKER - TEST SERVER")
    print("ws://localhost:8765")
    print("=" * 50)
    print()
    
    server = await asyncio.start_server(handle_client, 'localhost', 8765)
    
    async with server:
        print("[TEST] Server running. Open the HTML file in a browser.")
        print("[TEST] Press Ctrl+C to stop")
        print()
        await server.serve_forever()

if __name__ == '__main__':
    asyncio.run(main())