/**
 * interpolate.js — pill parameter interpolation + Q shaping + compile.
 *
 * pill: one template entry from heritage_designer_sections.json
 *   .sections: array of 6 { index, type, low_freq, low_gain, high_freq, high_gain }
 *
 * morph01: 0..1 (morph position)
 * qOffset01: -1..1 (Q offset, 0 = neutral)
 *
 * applyLiveQ shifts freq: at q=+1 → +50%, at q=-1 → -50%.
 * compileSections → number[6][5] via compileStage (shift=0).
 */

import { compileStage } from './emu_compile.js';

const FREQ_MIN = 0;
const FREQ_MAX = 127;

/**
 * Linear interpolation between low_* and high_* fields.
 * Returns array of 6 SectionParam objects:
 *   { index, type, freq, gain }
 */
export function sectionsAtMorph(pill, morph01) {
    const m = Math.max(0, Math.min(1, morph01));
    return pill.sections.slice(0, 6).map(sec => {
        const freq = sec.low_freq + m * (sec.high_freq - sec.low_freq);
        const gain = sec.low_gain + m * (sec.high_gain - sec.low_gain);
        return {
            index: sec.index,
            type:  sec.type,
            freq:  Math.round(freq),
            gain:  Math.round(gain),
        };
    });
}

/**
 * Apply live Q shaping: ±50% freq shift.
 * q=0 → no-op; q=+1 → freq * 1.5; q=-1 → freq * 0.5.
 * Clamps to [FREQ_MIN, FREQ_MAX].
 */
export function applyLiveQ(sectionParams, qOffset01) {
    const q = Math.max(-1, Math.min(1, qOffset01));
    if (q === 0) return sectionParams.map(s => ({ ...s }));
    return sectionParams.map(sec => {
        const shifted = sec.freq * (1.0 + 0.5 * q);
        const clamped = Math.round(Math.max(FREQ_MIN, Math.min(FREQ_MAX, shifted)));
        return { ...sec, freq: clamped };
    });
}

/**
 * Compile 6 section params to coefficients.
 * Returns number[6][5].
 * shift = 0 (heritage_designer paths use shift=0).
 */
export function compileSections(sectionParams, sr) {
    return sectionParams.map(sec =>
        compileStage(sec.type, sec.freq, sec.gain, 0, sr)
    );
}
