#!/usr/bin/env python3
"""
TRENCH Body Compiler Server
WebSocket bridge: receives body JSON from Body Baker, compiles to cartridge bytes
Targets: ws://localhost:8765
"""

import asyncio
import json
import struct
import math
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

SR = 39062.5


@dataclass
class StageState:
    kind: str
    pole_freq_hz: float
    pole_r: float
    zero_freq_hz: float
    zero_r: float


@dataclass
class CornerFrame:
    stages: List[StageState]


class MinifloatCodec:
    """16-bit minifloat encoder matching E-mu hardware spec"""
    
    @staticmethod
    def encode(value: float) -> int:
        """Encode float to 16-bit minifloat (4-bit exp, 12-bit mantissa)"""
        if value <= 0:
            return 0
        
        # Find exponent
        exp = math.floor(math.log2(value))
        if exp < -15:
            exp = -15
        if exp > 15:
            exp = 15
        
        # Compute mantissa (normalized to [1, 2))
        mantissa = value / (2 ** exp)
        mantissa_int = round((mantissa - 1) * 4096)
        
        if mantissa_int > 4095:
            mantissa_int = 4095
        
        biased_exp = exp + 15
        return (biased_exp << 12) | mantissa_int
    
    @staticmethod
    def decode(code: int) -> float:
        """Decode 16-bit minifloat to float"""
        if code == 0:
            return 0.0
        if code == 0xFFFF:
            return 1e10  # Max constant sentinel
        
        exp = (code >> 12) - 15
        mantissa = (code & 0xFFF) / 4096 + 1
        return mantissa * (2 ** exp)


class KernelForm:
    """Dave Rossum DF2T kernel form coefficients"""
    
    K = 2.0  # Identity constant
    
    @staticmethod
    def from_pole_zero(pole_freq: float, pole_r: float, zero_freq: float, zero_r: float) -> Tuple[float, float, float, float, float]:
        """Convert pole/zero to kernel form [c0, c1, c2, c3, c4]"""
        # These formulas derived from TRENCH RE codebase
        theta = 2 * math.pi * pole_freq / SR
        r = pole_r
        
        # Numerator (zeros)
        a = zero_r * math.cos(theta)
        b = zero_r * math.sin(theta)
        
        zr = zero_r
        zr2 = zr * zr
        
        b0 = 1 - 2 * a * zr + zr2
        b1 = 2 * a - 2 * zr
        b2 = 1
        
        # Denominator (poles)  
        a1 = -2 * r * math.cos(theta)
        a2 = r * r
        
        c0 = b0
        c1 = b1 + a1
        c2 = b2 + a1 * b1 + a2 * b0
        c3 = a1 * b2 + a2 * b1
        c4 = a2 * b2
        
        return c0, c1, c2, c3, c4


class BodyCompiler:
    """Compiles body JSON to compiled-v1 cartridge format"""
    
    CORNERS = ['M0_Q0', 'M0_Q100', 'M100_Q0', 'M100_Q100']
    
    def __init__(self):
        self.minifloat = MinifloatCodec()
        self.kernel = KernelForm()
    
    def compile(self, body: dict) -> dict:
        """Compile body JSON to cartridge format"""
        result = {
            'format': 'compiled-v1',
            'name': body.get('name', 'UNKNOWN'),
            'sample_rate': SR,
            'stage_count': body.get('stage_count', 6),
            'corners': {}
        }
        
        for corner_name in self.CORNERS:
            if corner_name not in body.get('corners', {}):
                continue
            
            corner_data = body['corners'][corner_name]
            compiled_stages = []
            
            for i, stage in enumerate(corner_data.get('stages', [])):
                if stage.get('kind') == 'off' or stage.get('pole_freq_hz', 0) == 0:
                    # Bypass stage: all coefficients = 0 (unity gain)
                    compiled = [0] * 5
                else:
                    # Compute kernel form coefficients
                    c0, c1, c2, c3, c4 = self.kernel.from_pole_zero(
                        stage['pole_freq_hz'],
                        stage['pole_r'],
                        stage['zero_freq_hz'],
                        stage['zero_r']
                    )
                    
                    # Minifloat encode
                    compiled = [
                        self.minifloat.encode(c0),
                        self.minifloat.encode(c1),
                        self.minifloat.encode(c2),
                        self.minifloat.encode(c3),
                        self.minifloat.encode(c4)
                    ]
                
                compiled_stages.append(compiled)
            
            result['corners'][corner_name] = compiled_stages
        
        return result
    
    def to_bytes(self, compiled: dict) -> bytes:
        """Serialize compiled body to binary cartridge format"""
        data = bytearray()
        
        # Header
        name = compiled['name'].encode('ascii')[:32].ljust(32, b'\x00')
        data.extend(name)
        
        data.extend(struct.pack('<I', 0x0100))  # Version
        data.extend(struct.pack('<d', compiled['sample_rate']))
        data.extend(struct.pack('<H', compiled['stage_count']))
        
        # Corner banks (4 corners × 6 stages × 5 uint16)
        for corner_name in self.CORNERS:
            if corner_name in compiled['corners']:
                stages = compiled['corners'][corner_name]
                for stage_coeffs in stages:
                    for coef in stage_coeffs:
                        data.extend(struct.pack('<H', coef))
            else:
                # Missing corner = zeros
                data.extend(b'\x00' * (6 * 5 * 2))
        
        return bytes(data)


class TrenchServer:
    """WebSocket server for TRENCH Body Baker"""
    
    def __init__(self, host: str = 'localhost', port: int = 8765):
        self.host = host
        self.port = port
        self.compiler = BodyCompiler()
        self.clients: List = []
    
    async def handle_client(self, reader: asyncio.StreamReader, writer: asyncio.StreamWriter):
        """Handle a single WebSocket client connection"""
        addr = writer.get_extra_info('peername')
        print(f"[WS] Client connected: {addr}")
        self.clients.append(writer)
        
        buffer = b''
        
        try:
            while True:
                data = await reader.read(65536)
                if not data:
                    break
                
                buffer += data
                
                # Process complete WebSocket frames
                while len(buffer) >= 2:
                    # Parse frame header
                    first_byte = buffer[0]
                    second_byte = buffer[1]
                    
                    opcode = first_byte & 0x0F
                    masked = (second_byte & 0x80) != 0
                    
                    payload_len = second_byte & 0x7F
                    header_len = 2
                    
                    if masked:
                        header_len += 4
                    
                    if payload_len == 126:
                        if len(buffer) < 4:
                            break
                        payload_len = struct.unpack('>H', buffer[2:4])[0]
                        header_len = 4
                    elif payload_len == 127:
                        if len(buffer) < 10:
                            break
                        payload_len = struct.unpack('>Q', buffer[2:10])[0]
                        header_len = 10
                    
                    if len(buffer) < header_len + payload_len:
                        break
                    
                    payload = buffer[header_len:header_len + payload_len]
                    buffer = buffer[header_len + payload_len:]
                    
                    # Handle close
                    if opcode == 0x8:
                        print(f"[WS] Client disconnected: {addr}")
                        break
                    
                    # Handle text (body JSON)
                    if opcode == 0x1:
                        try:
                            body = json.loads(payload.decode('utf-8'))
                            print(f"[WS] Received body: {body.get('name', 'unnamed')}")
                            
                            # Compile the body
                            compiled = self.compiler.compile(body)
                            cartridge_bytes = self.compiler.to_bytes(compiled)
                            
                            print(f"[WS] Compiled: {len(cartridge_bytes)} bytes")
                            
                            # Echo back as confirmation
                            response = json.dumps({
                                'status': 'compiled',
                                'cartridge_bytes': len(cartridge_bytes),
                                'corners': list(compiled['corners'].keys())
                            }).encode('utf-8')
                            
                            # Send WebSocket response
                            response_frame = self._frame_response(response)
                            writer.write(response_frame)
                            await writer.drain()
                            
                        except json.JSONDecodeError as e:
                            print(f"[WS] JSON parse error: {e}")
                    
                    # Handle pong
                    elif opcode == 0xA:
                        pass
        
        except Exception as e:
            print(f"[WS] Error: {e}")
        
        finally:
            if writer in self.clients:
                self.clients.remove(writer)
            writer.close()
            await writer.wait_closed()
            print(f"[WS] Client closed: {addr}")
    
    def _frame_response(self, payload: bytes) -> bytes:
        """Create WebSocket frame for response"""
        frame = bytearray()
        frame.append(0x81)  # FIN + text opcode
        
        if len(payload) < 126:
            frame.append(len(payload))
        elif len(payload) < 65536:
            frame.extend([0x7E, *struct.pack('>H', len(payload))])
        else:
            frame.extend([0x7F, *struct.pack('>Q', len(payload))])
        
        return bytes(frame) + payload
    
    async def start(self):
        """Start the WebSocket server"""
        server = await asyncio.start_server(
            self.handle_client,
            self.host,
            self.port
        )
        
        addr = server.sockets[0].getsockname()
        print(f"[WS] TRENCH Body Compiler listening on ws://{addr[0]}:{addr[1]}")
        print("[WS] Ready to receive body JSON from Body Baker")
        
        async with server:
            await server.serve_forever()


def main():
    print("=" * 60)
    print("TRENCH BODY COMPILER SERVER")
    print("WebSocket bridge for Body Baker → TRENCH engine")
    print("=" * 60)
    print()
    
    server = TrenchServer(host='localhost', port=8765)
    
    try:
        asyncio.run(server.start())
    except KeyboardInterrupt:
        print("\n[WS] Shutdown requested")


if __name__ == '__main__':
    main()