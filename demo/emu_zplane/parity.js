/**
 * parity.js — Cross-language parity check against hedz_golden.json.
 * Tolerance: maxAbsErr < 1e-4
 *
 * NOTE: The hedz baked ROM is a MorphDesigner-only template — Q is not
 * differentiated in the compiled kernel. hedz_golden.json therefore has
 * degenerate Q corners: M0_Q0 === M0_Q100 and M100_Q0 === M100_Q100.
 * This is confirmed in trench-core/src/hedz_golden.rs (build_cartridge
 * returns [m0, m1, m0, m1]).  Only the 2 morph endpoints at Q=0 are
 * tested here; the live-Q path is exercised separately by the slider in
 * the demo UI (applyLiveQ runs on top of the static kernel).
 */

import { renderOffline, makeImpulse } from './render_audio.js';
import { Cascade } from './cascade.js';

const BOOST = 4.0;

/**
 * Run the 4-corner parity check.
 * Returns array of { label, maxAbsErr, ok } for each corner.
 *
 * @param {object} interpolateMod  — module with sectionsAtMorph, applyLiveQ, compileSections
 * @param {number} sr              — sample rate (should be 44100 for golden)
 */
export async function runParity(interpolateMod, sr) {
    const { sectionsAtMorph, applyLiveQ, compileSections } = interpolateMod;

    // Fetch data
    const [goldenResp, sectionsResp] = await Promise.all([
        fetch('./data/hedz_golden.json'),
        fetch('./data/heritage_designer_sections.json'),
    ]);
    if (!goldenResp.ok)
        throw new Error(`HTTP ${goldenResp.status} fetching hedz_golden.json`);
    if (!sectionsResp.ok)
        throw new Error(`HTTP ${sectionsResp.status} fetching heritage_designer_sections.json`);
    const golden   = await goldenResp.json();
    const sections = await sectionsResp.json();

    // Find hedz pill
    const hedz = sections.templates.find(t => t.name === 'hedz');
    if (!hedz) throw new Error('hedz template not found in heritage_designer_sections.json');

    // Only 2 corners: M0 and M100 at Q=0.
    // Q corners are degenerate in the hedz baked ROM — see file header.
    const corners = [
        { label: 'M0_Q0',   morph: 0, q: 0 },
        { label: 'M100_Q0', morph: 1, q: 0 },
    ];

    const results = [];
    for (const { label, morph, q } of corners) {
        // Build target kernels
        const sectionParams = sectionsAtMorph(hedz, morph);
        const qed           = applyLiveQ(sectionParams, q);
        const target        = compileSections(qed, sr);

        // Build cascade (fresh per corner for parity — static coefficients, no ramp)
        const cascade = new Cascade();
        // rampSamples=1 → deltas are exactly (target - current) → snaps in 1 sample
        cascade.setTargets(target, 1);
        cascade.setBoost(BOOST, 1);

        // Render impulse
        const stimulus = makeImpulse(golden.len, sr);
        const out      = renderOffline(cascade, stimulus);

        // Compare
        const ref = golden[label];
        let maxAbsErr = 0;
        for (let i = 0; i < golden.len; i++) {
            const err = Math.abs(out[i] - ref[i]);
            if (err > maxAbsErr) maxAbsErr = err;
        }

        results.push({
            label,
            maxAbsErr,
            ok: maxAbsErr < 1e-4,
        });
    }
    return results;
}
