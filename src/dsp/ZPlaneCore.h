#pragma once
// =============================================================================
// ZPlaneCore.h — 14-Pole Z-Plane Filter Engine with Pole Morphing
// =============================================================================
// Core DSP engine for FIELD Z-Plane morphing filter.
//
// Key responsibilities:
//   - 7 biquad stages (14 poles total) in parallel or cascade topology
//   - Pole-parameter interpolation as primary morph domain
//   - Dual-path crossfade when topology/type differs between frames
//   - Per-stage bandwidth loudness compensation folded into gain
//   - Brute-force 7! pole assignment to minimize perceptual distance
//   - Saturation in feedback path (Dattorro soft clipper)
//   - Control-rate coefficient rebuild (every 16 samples)
//
// Architecture patterns:
//   - Audio thread safe: no allocations, no locks, no logging
//   - Double precision for coefficients/state, float for audio buffers
//   - Birth/death handling for pole activation transitions
//   - Inactive parallel stages output 0 (not bypass)
//
// Integration:
//   - Used by PluginProcessor for L/R channel processing
//   - Frames created by ZPlaneDesigner factory functions
//   - See docs/dsp-spec.md for complete mathematical reference
//
// Key sections:
//   - Pole struct (32-47): Semantic pole parameters (the morph domain)
//   - Biquad struct (53-66): Runtime coefficients with effectiveGain
//   - FilterFrame struct (79-89): Snapshot of pole positions + topology
//   - ZPlaneCore class (95-472): Main DSP engine
//     - Configuration (101-109): Sample rate and gain staging
//     - Frame setup (115-128): Set A/B frames, trigger assignment
//     - Morph control (133-137): Set interpolation position
//     - Processing (143-178): Sample/block processing with dual-path logic
//     - Pole interpolation (204-238): Birth/death/normal morph cases
//     - BW compensation (244-250): Per-stage loudness compensation
//     - Coefficient rebuild (256-288): Control-rate updates
//     - Pole assignment (351-374): Brute-force 7! permutation search
//     - Signal path (396-439): Biquad processing, parallel/cascade topology
//
// Usage example in PluginProcessor.cpp:62-70, 94-97, 155-156
// =============================================================================

#include <array>
#include <cmath>
#include <algorithm>

namespace field::dsp {

// =============================================================================
// Constants
// =============================================================================

static constexpr double PI = 3.14159265358979323846;
static constexpr int NUM_STAGES = 7;
static constexpr int CONTROL_RATE = 16;  // Rebuild coeffs every N samples
static constexpr double RADIUS_EPSILON = 1e-6;
static constexpr double BW_REF = 300.0;  // Reference bandwidth for loudness comp
static constexpr double BW_COMP_MAX = 12.0;  // Clamp compensation gain

// =============================================================================
// Pole — Semantic pole parameters (the morph domain)
// =============================================================================

struct Pole {
    double freq = 1000.0;    // Hz (20 - 20000)
    double radius = 0.0;     // 0.0 - 0.999 (0 = inactive/muted)
    double weight = 1.0;     // Per-stage amplitude (parallel mode)
    int id = -1;             // Identity for tracking through morphs

    bool isActive() const { return radius > RADIUS_EPSILON; }

    double distanceTo(const Pole& other) const {
        if (!isActive() && !other.isActive()) return 0.0;
        if (!isActive() || !other.isActive()) return 1e6;
        double d_freq = std::abs(std::log(freq) - std::log(other.freq));
        double d_radius = std::abs(std::log(radius) - std::log(other.radius));
        return d_freq + d_radius;
    }
};

// =============================================================================
// Biquad — Runtime coefficients
// =============================================================================

struct Biquad {
    double b0 = 1.0, b1 = 0.0, b2 = 0.0;
    double a1 = 0.0, a2 = 0.0;
    double effectiveGain = 1.0;  // weight * bwComp folded together
    double radius = 0.0;
    bool active = false;

    static Biquad muted() {
        Biquad b;
        b.active = false;
        b.effectiveGain = 0.0;
        return b;
    }
};

// =============================================================================
// Topology & FilterType
// =============================================================================

enum class Topology { PARALLEL, CASCADE };
enum class FilterType { RESONATOR, ALLPASS, LOWPASS };

// =============================================================================
// FilterFrame — A snapshot of pole positions + topology
// =============================================================================

struct FilterFrame {
    std::array<Pole, NUM_STAGES> poles{};
    Topology topology = Topology::PARALLEL;
    FilterType filterType = FilterType::RESONATOR;

    int activeCount() const {
        int count = 0;
        for (const auto& p : poles) if (p.isActive()) count++;
        return count;
    }
};

// =============================================================================
// ZPlaneCore — The DSP Engine with Dual-Path Crossfade
// =============================================================================

class ZPlaneCore {
public:
    // =========================================================================
    // Configuration
    // =========================================================================

    void setSampleRate(double sr) {
        sampleRate_ = sr;
        nyquist_ = sr * 0.45;
    }

    void setPreGain(float g) { preGain_ = g; }
    void setDrive(float d) { drive_ = d; }
    void setPostGain(float g) { postGain_ = g; }
    void setTrim(float t) { trim_ = t; }

    // =========================================================================
    // Frame Setup
    // =========================================================================

    void setFrameA(const FilterFrame& frame) {
        frameA_ = frame;
        computePoleAssignment();
        needsDualPath_ = checkNeedsDualPath();
        rebuildAllCoefficients();
    }

    void setFrameB(const FilterFrame& frame) {
        frameB_ = frame;
        computePoleAssignment();
        needsDualPath_ = checkNeedsDualPath();
        rebuildAllCoefficients();
    }

    // =========================================================================
    // Morph Control
    // =========================================================================

    void setMorph(float t) {
        morphT_ = std::clamp(static_cast<double>(t), 0.0, 1.0);
    }

    float getMorph() const { return static_cast<float>(morphT_); }

    // =========================================================================
    // Audio Processing
    // =========================================================================

    float processSample(float input) {
        if (++sampleCounter_ >= CONTROL_RATE) {
            sampleCounter_ = 0;
            rebuildAllCoefficients();
        }

        double x = input * preGain_ * drive_;
        double y;

        if (needsDualPath_) {
            // True A vs B crossfade: run both paths with their own coefficients
            double yA = processPath(x, stagesA_, stateA_, frameA_.topology);
            double yB = processPath(x, stagesB_, stateB_, frameB_.topology);
            y = (1.0 - morphT_) * yA + morphT_ * yB;
        } else {
            // Same topology/type: single interpolated path
            y = processPath(x, stagesA_, stateA_, frameA_.topology);
        }

        y = softclip(y * postGain_) * trim_;
        return static_cast<float>(y);
    }

    void processBlock(float* buffer, int numSamples) {
        for (int i = 0; i < numSamples; ++i) {
            buffer[i] = processSample(buffer[i]);
        }
    }

    void reset() {
        for (int i = 0; i < NUM_STAGES; ++i) {
            stateA_[i] = {};
            stateB_[i] = {};
        }
        sampleCounter_ = 0;
    }

private:
    // =========================================================================
    // Biquad State
    // =========================================================================

    struct BiquadState {
        double x1 = 0, x2 = 0;
        double y1 = 0, y2 = 0;
        double y1_prev = 0;
    };

    // =========================================================================
    // Check if dual path needed
    // =========================================================================

    bool checkNeedsDualPath() const {
        return (frameA_.topology != frameB_.topology) ||
               (frameA_.filterType != frameB_.filterType);
    }

    // =========================================================================
    // Pole Interpolation
    // =========================================================================

    Pole interpolatePole(const Pole& a, const Pole& b, double t) const {
        Pole result;

        if (!a.isActive() && !b.isActive()) {
            result.radius = 0.0;
            result.freq = 1000.0;
            result.weight = 0.0;
            return result;
        }

        // Birth
        if (!a.isActive() && b.isActive()) {
            result.freq = b.freq;
            result.radius = std::pow(b.radius, t);
            result.weight = t * b.weight;
            return result;
        }

        // Death
        if (a.isActive() && !b.isActive()) {
            result.freq = a.freq;
            result.radius = std::pow(a.radius, 1.0 - t);
            result.weight = (1.0 - t) * a.weight;
            return result;
        }

        // Normal morph
        result.freq = a.freq * std::pow(b.freq / a.freq, t);
        result.freq = std::clamp(result.freq, 20.0, nyquist_);
        result.radius = std::pow(a.radius, 1.0 - t) * std::pow(b.radius, t);
        result.radius = std::clamp(result.radius, 0.0, 1.0 - RADIUS_EPSILON);
        result.weight = (1.0 - t) * a.weight + t * b.weight;

        return result;
    }

    // =========================================================================
    // Per-stage bandwidth compensation
    // =========================================================================

    double computeStageBwGain(double radius) const {
        if (radius <= RADIUS_EPSILON) return 1.0;
        double bwHz = -std::log(radius) * sampleRate_ / PI;
        if (bwHz < 1.0) bwHz = 1.0;
        double g = std::sqrt(BW_REF / bwHz);
        return std::clamp(g, 1.0, BW_COMP_MAX);
    }

    // =========================================================================
    // Coefficient Rebuilding
    // =========================================================================

    void rebuildAllCoefficients() {
        if (needsDualPath_) {
            // Build separate coefficient sets for A and B
            rebuildCoefficientsForFrame(frameA_, stagesA_, 0.0);
            rebuildCoefficientsForFrame(frameB_, stagesB_, 1.0);
        } else {
            // Single interpolated set
            rebuildInterpolatedCoefficients();
        }
    }

    void rebuildCoefficientsForFrame(const FilterFrame& frame,
                                      std::array<Biquad, NUM_STAGES>& stages,
                                      double t) {
        for (int i = 0; i < NUM_STAGES; ++i) {
            const Pole& pole = frame.poles[i];
            stages[i] = poleToCoefficients(pole, frame.filterType);
        }
    }

    void rebuildInterpolatedCoefficients() {
        for (int i = 0; i < NUM_STAGES; ++i) {
            int mappedIdx = assignment_[i];
            Pole poleA = frameA_.poles[i];
            Pole poleB = (mappedIdx >= 0) ? frameB_.poles[mappedIdx] : Pole{};

            if (birthMask_[i]) poleA.radius = 0.0;
            if (deathMask_[i]) poleB.radius = 0.0;

            Pole interpolated = interpolatePole(poleA, poleB, morphT_);
            stagesA_[i] = poleToCoefficients(interpolated, frameA_.filterType);
        }
    }

    // =========================================================================
    // Pole to Coefficients (with per-stage BW comp)
    // =========================================================================

    Biquad poleToCoefficients(const Pole& pole, FilterType type) const {
        Biquad b;

        if (!pole.isActive()) {
            b.active = false;
            b.b0 = 0.0; b.b1 = 0.0; b.b2 = 0.0;
            b.a1 = 0.0; b.a2 = 0.0;
            b.effectiveGain = 0.0;
            return b;
        }

        double theta = 2.0 * PI * pole.freq / sampleRate_;
        double R = pole.radius;

        b.a1 = -2.0 * R * std::cos(theta);
        b.a2 = R * R;
        b.radius = R;
        b.active = true;

        // Per-stage BW compensation folded into effective gain
        double bwGain = computeStageBwGain(R);
        b.effectiveGain = pole.weight * bwGain;

        switch (type) {
            case FilterType::RESONATOR: {
                double gain = 1.0 - R;
                b.b0 = gain;
                b.b1 = 0.0;
                b.b2 = -gain;
                break;
            }
            case FilterType::ALLPASS: {
                b.b0 = R * R;
                b.b1 = -2.0 * R * std::cos(theta);
                b.b2 = 1.0;
                break;
            }
            case FilterType::LOWPASS: {
                double K = std::tan(PI * pole.freq / sampleRate_);
                double Q = 1.0 / (2.0 * (1.0 - R));
                double norm = 1.0 / (1.0 + K / Q + K * K);
                b.b0 = K * K * norm;
                b.b1 = 2.0 * b.b0;
                b.b2 = b.b0;
                b.a1 = 2.0 * (K * K - 1.0) * norm;
                b.a2 = (1.0 - K / Q + K * K) * norm;
                break;
            }
        }

        return b;
    }

    // =========================================================================
    // Brute-force 7! Pole Assignment
    // =========================================================================

    void computePoleAssignment() {
        // Generate all permutations and find minimum cost
        int perm[NUM_STAGES] = {0, 1, 2, 3, 4, 5, 6};
        int bestPerm[NUM_STAGES];
        double bestCost = 1e18;

        do {
            double cost = 0.0;
            for (int i = 0; i < NUM_STAGES; ++i) {
                cost += frameA_.poles[i].distanceTo(frameB_.poles[perm[i]]);
            }
            if (cost < bestCost) {
                bestCost = cost;
                for (int i = 0; i < NUM_STAGES; ++i) bestPerm[i] = perm[i];
            }
        } while (std::next_permutation(perm, perm + NUM_STAGES));

        // Store best assignment
        for (int i = 0; i < NUM_STAGES; ++i) {
            assignment_[i] = bestPerm[i];
            birthMask_[i] = !frameA_.poles[i].isActive() && frameB_.poles[bestPerm[i]].isActive();
            deathMask_[i] = frameA_.poles[i].isActive() && !frameB_.poles[bestPerm[i]].isActive();
        }
    }

    // =========================================================================
    // Saturation
    // =========================================================================

    static double saturate(double x) {
        if (x > 1.0) return 0.666666666666667;
        if (x < -1.0) return -0.666666666666667;
        return x - (x * x * x) / 3.0;
    }

    static double softclip(double x) {
        if (x > 1.5) return 1.0;
        if (x < -1.5) return -1.0;
        return x - (x * x * x) / 6.75;
    }

    // =========================================================================
    // Process single biquad
    // =========================================================================

    double processBiquad(const Biquad& b, BiquadState& s, double x) {
        double y = b.b0 * x + b.b1 * s.x1 + b.b2 * s.x2
                 - b.a1 * s.y1 - b.a2 * s.y2;

        s.x2 = s.x1;
        s.x1 = x;
        s.y1_prev = s.y1;
        s.y1 = saturate(y);
        s.y2 = s.y1_prev;

        return y;
    }

    // =========================================================================
    // Process path (parallel or cascade)
    // =========================================================================

    double processPath(double x,
                       std::array<Biquad, NUM_STAGES>& stages,
                       std::array<BiquadState, NUM_STAGES>& state,
                       Topology topology) {
        if (topology == Topology::PARALLEL) {
            double sum = 0.0;
            double gainSumSq = 0.0;

            for (int i = 0; i < NUM_STAGES; ++i) {
                if (!stages[i].active) continue;
                double y = processBiquad(stages[i], state[i], x);
                double g = stages[i].effectiveGain;
                sum += y * g;
                gainSumSq += g * g;
            }

            double norm = std::sqrt(gainSumSq) + 1e-12;
            return sum / norm;
        } else {
            double y = x;
            for (int i = 0; i < NUM_STAGES; ++i) {
                if (!stages[i].active) continue;
                y = processBiquad(stages[i], state[i], y);
            }
            return y;
        }
    }

    // =========================================================================
    // State
    // =========================================================================

    FilterFrame frameA_{};
    FilterFrame frameB_{};
    double morphT_ = 0.0;
    bool needsDualPath_ = false;

    // Dual coefficient sets
    std::array<Biquad, NUM_STAGES> stagesA_{};
    std::array<Biquad, NUM_STAGES> stagesB_{};

    // Dual state sets
    std::array<BiquadState, NUM_STAGES> stateA_{};
    std::array<BiquadState, NUM_STAGES> stateB_{};

    // Pole assignment (brute-force result)
    int assignment_[NUM_STAGES] = {0, 1, 2, 3, 4, 5, 6};
    bool birthMask_[NUM_STAGES] = {false};
    bool deathMask_[NUM_STAGES] = {false};

    // Gain staging
    float preGain_ = 1.0f;
    float drive_ = 1.0f;
    float postGain_ = 4.0f;
    float trim_ = 0.7f;

    double sampleRate_ = 48000.0;
    double nyquist_ = 21600.0;
    int sampleCounter_ = 0;
};

} // namespace field::dsp
