/**
 * cascade.js — 12-stage serial DF2T biquad cascade.
 * Ported from trench-core/src/cascade.rs.
 *
 * DF2T equations (exact):
 *   y  = c0*x + w1
 *   w1 = c1*x - c3*y + w2
 *   w2 = c2*x - c4*y
 */

const TOTAL_STAGES = 12;
const NUM_STAGES   = 6;  // active
const PASSTHROUGH  = [1, 0, 0, 0, 0];

class BiquadState {
    constructor() {
        this.coeffs = [1, 0, 0, 0, 0];
        this.deltas = [0, 0, 0, 0, 0];
        this.w1 = 0.0;
        this.w2 = 0.0;
    }

    setTarget(target, rampSamples) {
        const n = Math.max(1, rampSamples);
        for (let i = 0; i < 5; i++) {
            this.deltas[i] = (target[i] - this.coeffs[i]) / n;
        }
    }

    /** Returns [output, unstable]. Ramps coefficients, then computes DF2T. */
    processSample(x) {
        // Ramp
        this.coeffs[0] += this.deltas[0];
        this.coeffs[1] += this.deltas[1];
        this.coeffs[2] += this.deltas[2];
        this.coeffs[3] += this.deltas[3];
        this.coeffs[4] += this.deltas[4];

        // DF2T
        let y = this.coeffs[0] * x + this.w1;
        if (!isFinite(y)) {
            y = 0.0;
            this.w1 = 0.0;
            this.w2 = 0.0;
            this.deltas = [0, 0, 0, 0, 0];
            return [y, true];
        }
        this.w1 = this.coeffs[1] * x - this.coeffs[3] * y + this.w2;
        this.w2 = this.coeffs[2] * x - this.coeffs[4] * y;
        const unstable = !isFinite(this.w1) || !isFinite(this.w2);
        if (unstable) {
            this.w1 = 0.0;
            this.w2 = 0.0;
            this.deltas = [0, 0, 0, 0, 0];
        }
        return [y, unstable];
    }
}

export class Cascade {
    constructor() {
        this._stages = [];
        for (let i = 0; i < TOTAL_STAGES; i++) {
            this._stages.push(new BiquadState());
        }
        this._boost      = 1.0;
        this._boostDelta = 0.0;
        this._instabilityFlag = false;
    }

    /**
     * Set target coefficients for 6 active stages.
     * target: number[6][5] — stages 0..5.
     * Stages 6..11 are set to passthrough [1,0,0,0,0].
     * rampSamples: number of samples to ramp over.
     */
    setTargets(target, rampSamples) {
        for (let i = 0; i < NUM_STAGES; i++) {
            this._stages[i].setTarget(target[i], rampSamples);
        }
        for (let i = NUM_STAGES; i < TOTAL_STAGES; i++) {
            this._stages[i].setTarget(PASSTHROUGH, rampSamples);
        }
    }

    /** Set boost target, ramped over rampSamples. */
    setBoost(target, rampSamples) {
        const n = Math.max(1, rampSamples);
        this._boostDelta = (target - this._boost) / n;
    }

    /** Process one sample through all 12 stages + boost. */
    processSample(x) {
        let v = x;
        for (const stage of this._stages) {
            const [next, unstable] = stage.processSample(v);
            if (unstable) {
                this._instabilityFlag = true;
                return 0.0;
            }
            v = next;
        }
        this._boost += this._boostDelta;
        if (!isFinite(this._boost)) {
            this._boost = 1.0;
            this._boostDelta = 0.0;
            this._instabilityFlag = true;
            return 0.0;
        }
        v *= this._boost;
        if (!isFinite(v)) {
            this._instabilityFlag = true;
            return 0.0;
        }
        return v;
    }

    /** Process Float32Array in-place, 32 samples at a time. */
    processBlockMono(samples) {
        for (let i = 0; i < samples.length; i++) {
            samples[i] = this.processSample(samples[i]);
        }
    }

    /** Returns and clears the instability flag. */
    takeInstabilityFlag() {
        const f = this._instabilityFlag;
        this._instabilityFlag = false;
        return f;
    }

    /** Reset all delay state. */
    reset() {
        for (const stage of this._stages) {
            stage.w1 = 0.0;
            stage.w2 = 0.0;
            stage.deltas = [0, 0, 0, 0, 0];
        }
        this._boostDelta = 0.0;
        this._instabilityFlag = false;
    }
}
