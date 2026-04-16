/**
 * plot.js — Canvas 2D response and impulse drawing.
 * Bassbox-303 palette: bg #111, grid #2a3, trace #fb3, axes #555, labels #8b8.
 */

const DB_MIN   = -60;
const DB_MAX   = 12;
const FREQ_MIN = 20;

/**
 * Draw frequency response on canvas.
 * magDb: Float32Array of FFT magnitude in dB (DC..Nyquist).
 * sr: sample rate (determines Nyquist cap).
 */
export function drawResponse(canvas, magDb, sr) {
    const ctx    = canvas.getContext('2d');
    const W      = canvas.width;
    const H      = canvas.height;
    const freqMax = Math.min(sr / 2, 20000);
    const padL   = 42;
    const padR   = 10;
    const padT   = 10;
    const padB   = 24;
    const plotW  = W - padL - padR;
    const plotH  = H - padT - padB;

    ctx.fillStyle = '#111';
    ctx.fillRect(0, 0, W, H);

    // Frequency axis: log scale from FREQ_MIN to freqMax
    const logMin = Math.log10(FREQ_MIN);
    const logMax = Math.log10(freqMax);
    const logRange = logMax - logMin;

    function freqToX(f) {
        return padL + ((Math.log10(f) - logMin) / logRange) * plotW;
    }
    function dbToY(db) {
        const t = (DB_MAX - db) / (DB_MAX - DB_MIN);
        return padT + t * plotH;
    }

    // Grid lines — frequency
    ctx.strokeStyle = '#2a3';
    ctx.lineWidth   = 0.5;
    const freqMarks = [20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000];
    for (const f of freqMarks) {
        if (f > freqMax) break;
        const x = freqToX(f);
        ctx.beginPath();
        ctx.moveTo(x, padT);
        ctx.lineTo(x, padT + plotH);
        ctx.stroke();
    }

    // Grid lines — dB
    const dbMarks = [-60, -48, -36, -24, -12, 0, 12];
    for (const db of dbMarks) {
        const y = dbToY(db);
        ctx.beginPath();
        ctx.moveTo(padL, y);
        ctx.lineTo(padL + plotW, y);
        ctx.stroke();
    }

    // 0 dB line brighter
    ctx.strokeStyle = '#555';
    ctx.lineWidth = 1;
    const y0 = dbToY(0);
    ctx.beginPath();
    ctx.moveTo(padL, y0);
    ctx.lineTo(padL + plotW, y0);
    ctx.stroke();

    // Axes
    ctx.strokeStyle = '#555';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, padT);
    ctx.lineTo(padL, padT + plotH);
    ctx.lineTo(padL + plotW, padT + plotH);
    ctx.stroke();

    // Labels
    ctx.fillStyle = '#8b8';
    ctx.font = '9px monospace';
    ctx.textAlign = 'center';
    for (const f of freqMarks) {
        if (f > freqMax) break;
        const x   = freqToX(f);
        const lbl = f >= 1000 ? `${f / 1000}k` : `${f}`;
        ctx.fillText(lbl, x, padT + plotH + 14);
    }
    ctx.textAlign = 'right';
    for (const db of dbMarks) {
        ctx.fillText(`${db}`, padL - 4, dbToY(db) + 3);
    }

    // Trace
    // magDb is indexed 0..N/2, mapping to 0..sr/2
    const nyquist  = sr / 2;
    const nBins    = magDb.length;

    ctx.strokeStyle = '#fb3';
    ctx.lineWidth   = 1.5;
    ctx.beginPath();
    let started = false;
    for (let i = 1; i < nBins; i++) {
        const freq = (i / (nBins - 1)) * nyquist;
        if (freq < FREQ_MIN || freq > freqMax) continue;
        const x  = freqToX(freq);
        const db = Math.max(DB_MIN, Math.min(DB_MAX, magDb[i]));
        const y  = dbToY(db);
        if (!started) { ctx.moveTo(x, y); started = true; }
        else ctx.lineTo(x, y);
    }
    ctx.stroke();
}

/**
 * Draw first 256 samples of impulse response.
 */
export function drawImpulse(canvas, samples) {
    const ctx = canvas.getContext('2d');
    const W   = canvas.width;
    const H   = canvas.height;
    const n   = Math.min(256, samples.length);
    const padL = 6;
    const padR = 6;
    const padT = 6;
    const padB = 6;
    const plotW = W - padL - padR;
    const plotH = H - padT - padB;

    ctx.fillStyle = '#111';
    ctx.fillRect(0, 0, W, H);

    // Find range
    let mn = -1, mx = 1;
    for (let i = 0; i < n; i++) {
        if (samples[i] < mn) mn = samples[i];
        if (samples[i] > mx) mx = samples[i];
    }
    const range = mx - mn || 1;

    // Zero line
    ctx.strokeStyle = '#555';
    ctx.lineWidth   = 0.5;
    const y0 = padT + plotH * (1 - (0 - mn) / range);
    ctx.beginPath();
    ctx.moveTo(padL, y0);
    ctx.lineTo(padL + plotW, y0);
    ctx.stroke();

    // Trace
    ctx.strokeStyle = '#fb3';
    ctx.lineWidth   = 1;
    ctx.beginPath();
    for (let i = 0; i < n; i++) {
        const x = padL + (i / (n - 1)) * plotW;
        const y = padT + plotH * (1 - (samples[i] - mn) / range);
        if (i === 0) ctx.moveTo(x, y);
        else ctx.lineTo(x, y);
    }
    ctx.stroke();
}
