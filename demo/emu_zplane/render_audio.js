/**
 * render_audio.js — stimulus generators and offline rendering.
 */

const BLOCK = 32;

/** Unit impulse of given length. */
export function makeImpulse(len, sr) {
    const buf = new Float32Array(len);
    buf[0] = 1.0;
    return buf;
}

/** Sine burst at hz for len samples. */
export function makeSineBurst(hz, len, sr) {
    const buf = new Float32Array(len);
    const w   = 2 * Math.PI * hz / sr;
    for (let i = 0; i < len; i++) buf[i] = Math.sin(w * i);
    return buf;
}

/** Deterministic white noise via LCG. seed is optional (default 42). */
export function makeWhiteNoise(len, sr, seed) {
    const buf = new Float32Array(len);
    let s = (seed !== undefined ? seed : 42) >>> 0;
    for (let i = 0; i < len; i++) {
        // Park-Miller LCG
        s = Math.imul(s, 1664525) + 1013904223;
        buf[i] = ((s >>> 1) / 0x40000000) - 1.0; // range [-1, 1)
    }
    return buf;
}

/** Short impulse train (click). */
export function makeClick(len, sr) {
    const buf    = new Float32Array(len);
    const period = Math.round(sr / 1000); // ~1 ms click rate
    for (let i = 0; i < len; i += period) {
        if (i < len) buf[i] = 1.0;
    }
    return buf;
}

/**
 * Run cascade over stimulus. cascade must be pre-configured (setTargets called).
 * Returns Float32Array output.
 */
export function renderOffline(cascade, stimulus) {
    const out = new Float32Array(stimulus.length);
    for (let i = 0; i < stimulus.length; i++) out[i] = stimulus[i];

    for (let i = 0; i < out.length; i += BLOCK) {
        const chunk = out.subarray(i, Math.min(i + BLOCK, out.length));
        cascade.processBlockMono(chunk);
    }
    return out;
}

/**
 * Wrap a mono Float32Array into an AudioBuffer.
 */
export function toAudioBuffer(ctx, mono) {
    const buf = ctx.createBuffer(1, mono.length, ctx.sampleRate);
    buf.copyToChannel(mono, 0);
    return buf;
}
