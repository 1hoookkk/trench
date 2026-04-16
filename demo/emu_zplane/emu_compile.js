/**
 * emu_compile.js — E-mu Z-plane firmware compile math.
 * Ported from tools/bake_hedz_const.py lines 65–172.
 * Mirrors pyruntime/heritage_coeffs.py type 1/2/3 recipes.
 * Math.fround used at all intermediate f32 storage points.
 */

export const COMBINE_K = 4.0;

export const SR_FAMILY = { 44100: 0, 48000: 1, 96000: 2, 192000: 3 };

export const FW_SCALE = [220, 220, 200, 177];
export const FW_BASE  = [18, 18, 4, 1];

// ---------------------------------------------------------------------------
// u16 minifloat → f32
// Mirrors trench-core/src/minifloat.rs and pyruntime/heritage_coeffs.py:_mf_decode
// ---------------------------------------------------------------------------

export function mfDecode(u16) {
    const packed = u16 & 0xFFFF;
    if (packed === 0xFFFF) return Math.fround(1.0);
    if (packed === 0x0000) return Math.fround(0.0);
    const u = packed + 1;
    const exp = ((u >> 12) & 0xF) - 15;
    const mant = u & 0xFFF;
    if (exp < -14) {
        return Math.fround(mant / 134217728.0); // 2^-27
    }
    const scaled = (mant | 0x1000) / 8192.0; // 2^-13
    return Math.fround(scaled * Math.pow(2.0, exp));
}

// ---------------------------------------------------------------------------
// 5 uint16 firmware words → direct-form kernel (c0..c4)
// ---------------------------------------------------------------------------

export function fwWordsToKernel(words) {
    const [w0, w1, w2, w3, w4] = words;
    const d0 = mfDecode(w0);
    const d1 = mfDecode(w1);
    const d2 = mfDecode(w2);
    const d3 = mfDecode(w3);
    const d4 = mfDecode(w4);
    return [
        Math.fround(d0 * COMBINE_K + d1),
        Math.fround(d1),
        Math.fround(d2 * COMBINE_K + d3),
        Math.fround(d3),
        Math.fround(d4),
    ];
}

// ---------------------------------------------------------------------------
// Frequency value from packed freq + sample rate index
// ---------------------------------------------------------------------------

export function fwFreqValue(packedFreq, sr) {
    const idx = SR_FAMILY[sr] !== undefined ? SR_FAMILY[sr] : 0;
    return Math.floor((FW_SCALE[idx] * packedFreq) / 128) + FW_BASE[idx];
}

function fwRadius(freqValue) {
    return Math.floor((freqValue * 124) / 256) + 118;
}

function fwGainOffset(packedGain, globalShift) {
    const raw = Math.floor((packedGain - 64) / 2) + (globalShift || 0);
    return Math.max(-32, Math.min(31, raw));
}

function clamp255(v) {
    return Math.max(0, Math.min(255, v));
}

// ---------------------------------------------------------------------------
// Type 1 kernel
// ---------------------------------------------------------------------------

export function type1Kernel(freqPacked, gainPacked, shift, sr) {
    const fv  = fwFreqValue(freqPacked, sr);
    const rad = fwRadius(fv);
    const go  = fwGainOffset(gainPacked, shift);
    const w0  = (fv << 8) >>> 0;
    const w1  = (clamp255(rad + go) << 8) >>> 0;
    const w2  = (fv << 8) >>> 0;
    const w3  = (clamp255(rad - go) << 8) >>> 0;
    const w4  = 0xE000;
    return fwWordsToKernel([w0, w1, w2, w3, w4]);
}

// ---------------------------------------------------------------------------
// Type 2 kernel
// ---------------------------------------------------------------------------

export function type2Kernel(freqPacked, gainPacked, shift, sr) {
    const idx = SR_FAMILY[sr] !== undefined ? SR_FAMILY[sr] : 0;
    const fv  = fwFreqValue(freqPacked, sr);
    const rad = fwRadius(fv);
    const go  = fwGainOffset(gainPacked, shift);
    let w0, w1, w4_val;
    if (idx < 2) {
        w0    = (0xEC << 8) >>> 0;
        w1    = (0xFF << 8) >>> 0;
        w4_val = fv + 0xF5;
    } else {
        w0    = (0xE1 << 8) >>> 0;
        w1    = (0xF0 << 8) >>> 0;
        w4_val = fv;
    }
    const w2 = (fv << 8) >>> 0;
    const w3 = (clamp255(rad - go) << 8) >>> 0;
    const w4 = (w4_val << 8) >>> 0;
    return fwWordsToKernel([w0, w1, w2, w3, w4]);
}

// ---------------------------------------------------------------------------
// Type 3 kernel
// ---------------------------------------------------------------------------

function type3FreqCompression(freqValue, shift) {
    if (freqValue > 0xDB && shift < 0) {
        freqValue = (((freqValue - 220) * (shift + 32)) >> 5) + 220;
    }
    return freqValue;
}

export function type3Kernel(freqPacked, gainPacked, shift, sr) {
    const idx  = SR_FAMILY[sr] !== undefined ? SR_FAMILY[sr] : 0;
    const fv   = fwFreqValue(freqPacked, sr);
    const go   = fwGainOffset(gainPacked, shift);
    const fvc  = type3FreqCompression(fv, shift);
    const rad  = fwRadius(fvc);
    const base = FW_BASE[idx];
    const w0   = (base << 8) >>> 0;
    const w1   = ((Math.floor((base * 124) / 256) + 150) << 8) >>> 0;
    const w2   = (fvc << 8) >>> 0;
    const w3   = (clamp255(rad - go) << 8) >>> 0;
    let w4;
    if (idx < 2) {
        // c4_raw = (fv_compressed - 18) * (-12) + (-8192)
        const c4_raw = (fvc - 18) * (-12) + (-8192);
        w4 = (c4_raw & 0xFFFF) >>> 0;
    } else {
        w4 = 0xE000;
    }
    return fwWordsToKernel([w0, w1, w2, w3, w4]);
}

// ---------------------------------------------------------------------------
// compile one stage → [c0..c4]. Type 0 = passthrough.
// ---------------------------------------------------------------------------

export function compileStage(typeId, freqPacked, gainPacked, shift, sr) {
    if (typeId <= 0) return [1, 0, 0, 0, 0];
    if (typeId === 1) return type1Kernel(freqPacked, gainPacked, shift, sr);
    if (typeId === 2) return type2Kernel(freqPacked, gainPacked, shift, sr);
    if (typeId === 3) return type3Kernel(freqPacked, gainPacked, shift, sr);
    // fallback: type 1
    return type1Kernel(freqPacked, gainPacked, shift, sr);
}
