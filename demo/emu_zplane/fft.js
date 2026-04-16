/**
 * fft.js — Hand-written radix-2 real FFT (N=1024).
 */

const FFT_N = 1024;

/**
 * In-place Cooley-Tukey FFT. re and im are Float64Arrays of length N.
 * Input must be N = power of 2.
 */
function fftInPlace(re, im, n) {
    // Bit-reverse permutation
    let j = 0;
    for (let i = 1; i < n; i++) {
        let bit = n >> 1;
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        if (i < j) {
            let t = re[i]; re[i] = re[j]; re[j] = t;
            t = im[i]; im[i] = im[j]; im[j] = t;
        }
    }
    // FFT butterfly
    for (let len = 2; len <= n; len <<= 1) {
        const half = len >> 1;
        const ang  = -2.0 * Math.PI / len;
        const wRe  = Math.cos(ang);
        const wIm  = Math.sin(ang);
        for (let i = 0; i < n; i += len) {
            let urRe = 1.0, urIm = 0.0;
            for (let k = 0; k < half; k++) {
                const a = i + k;
                const b = a + half;
                const tRe = urRe * re[b] - urIm * im[b];
                const tIm = urRe * im[b] + urIm * re[b];
                re[b] = re[a] - tRe;
                im[b] = im[a] - tIm;
                re[a] += tRe;
                im[a] += tIm;
                const newUrRe = urRe * wRe - urIm * wIm;
                urIm = urRe * wIm + urIm * wRe;
                urRe = newUrRe;
            }
        }
    }
}

/**
 * Real FFT of Float32Array x of length 1024.
 * Returns { re: Float64Array(513), im: Float64Array(513) }
 * (DC to Nyquist, inclusive).
 */
export function realFft(x) {
    const n   = FFT_N;
    const re  = new Float64Array(n);
    const im  = new Float64Array(n);
    for (let i = 0; i < n; i++) { re[i] = x[i]; }

    fftInPlace(re, im, n);

    // Return DC..Nyquist (n/2 + 1 bins)
    const half = n / 2 + 1;
    const outRe = new Float64Array(half);
    const outIm = new Float64Array(half);
    for (let i = 0; i < half; i++) {
        outRe[i] = re[i];
        outIm[i] = im[i];
    }
    return { re: outRe, im: outIm };
}

const DB_FLOOR = -80.0;

/**
 * Compute magnitude in dB from re/im arrays.
 * Returns Float32Array of same length. Floor: -80 dB.
 */
export function magnitudeDb(re, im) {
    const n  = re.length;
    const db = new Float32Array(n);
    for (let i = 0; i < n; i++) {
        const mag = Math.sqrt(re[i] * re[i] + im[i] * im[i]);
        db[i] = mag > 0 ? Math.max(DB_FLOOR, 20 * Math.log10(mag)) : DB_FLOOR;
    }
    return db;
}
