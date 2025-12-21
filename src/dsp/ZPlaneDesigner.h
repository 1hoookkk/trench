#pragma once
// =============================================================================
// ZPlaneDesigner.h — Pole and Frame Factory
// =============================================================================
// Factory functions for creating Pole objects and FilterFrames from
// frequency/bandwidth specifications.
//
// Key responsibilities:
//   - Convert semantic parameters (freq, bandwidth, weight) to Pole objects
//   - Provide bandwidth/radius/Q conversion utilities
//   - Define vowel formant presets (Peterson & Barney 1952 reference)
//   - Build complete FilterFrames for different filter types
//
// Frame builders:
//   - buildVowelFrame (84-106): PARALLEL RESONATOR for formant synthesis
//   - buildPhaserFrame (109-130): CASCADE ALLPASS for phase effects
//   - buildAcidFrame (133-157): CASCADE LOWPASS for resonant filters
//   - buildCustomFrame (160-179): Generic builder from pole array
//
// Conversions:
//   - bandwidthToRadius (40-43): BW (Hz) -> radius via R = exp(-PI*BW/sr)
//   - radiusToBandwidth (45-48): Inverse conversion
//   - qToBandwidth / bandwidthToQ (54-60): Q-factor conversions
//
// Vowel formants:
//   - Peterson & Barney 1952 formant frequencies
//   - Tight bandwidth (15-20 Hz) for E-mu "ring" character
//   - Per-formant weights for vowel character
//
// Integration:
//   - Creates FilterFrames consumed by ZPlaneCore
//   - Used by PluginProcessor.cpp:78-92 for FIELD parameter mapping
//   - Works in pole-parameter domain (not coefficients)
//
// See also:
//   - ZPlaneCore.h: Filter engine that consumes these frames
//   - docs/dsp-spec.md:63-118: Coefficient formulas and filter types
// =============================================================================

#include "ZPlaneCore.h"
#include <cmath>

namespace field::dsp {

class ZPlaneDesigner {
public:
    // =========================================================================
    // Pole Creation — The primary interface
    // =========================================================================

    // Create a pole from frequency and bandwidth (Hz)
    static Pole makePole(double freq, double bandwidth, double sr, double weight = 1.0) {
        Pole p;
        p.freq = std::clamp(freq, 20.0, sr * 0.45);
        p.radius = bandwidthToRadius(bandwidth, sr);
        p.weight = weight;
        return p;
    }

    // Create an inactive (muted) pole
    static Pole mutedPole() {
        Pole p;
        p.radius = 0.0;
        return p;
    }

    // =========================================================================
    // Bandwidth ↔ Radius Conversion
    // =========================================================================

    static double bandwidthToRadius(double bandwidth, double sr) {
        double R = std::exp(-PI * bandwidth / sr);
        return std::clamp(R, 0.0, 1.0 - RADIUS_EPSILON);
    }

    static double radiusToBandwidth(double radius, double sr) {
        if (radius <= 0.0) return sr;  // Infinite bandwidth
        return -std::log(radius) * sr / PI;
    }

    // =========================================================================
    // Q ↔ Bandwidth Conversion (for user-facing controls)
    // =========================================================================

    static double qToBandwidth(double Q, double freq) {
        return freq / Q;
    }

    static double bandwidthToQ(double bandwidth, double freq) {
        return freq / bandwidth;
    }
};

// =============================================================================
// VowelFormants — Peterson & Barney 1952
// =============================================================================

struct VowelFormants {
    double f1, f2, f3, f4;   // Frequencies (Hz)
    double bw;                // Bandwidth (Hz) - tight for E-mu ring
    double w1 = 1.0, w2 = 0.8, w3 = 0.5, w4 = 0.3;  // Per-formant weights

    static VowelFormants Ah() { return {730, 1090, 2440, 3500, 18, 1.0, 0.8, 0.5, 0.3}; }
    static VowelFormants Ee() { return {270, 2290, 3010, 3700, 16, 1.0, 0.9, 0.6, 0.3}; }
    static VowelFormants Oo() { return {300, 870, 2240, 3500, 15, 1.0, 0.8, 0.5, 0.3}; }
    static VowelFormants Eh() { return {530, 1840, 2480, 3500, 18, 1.0, 0.8, 0.5, 0.3}; }
    static VowelFormants Uh() { return {640, 1190, 2390, 3500, 20, 1.0, 0.8, 0.5, 0.3}; }
};

// =============================================================================
// Frame Builders — Create complete FilterFrames from specs
// =============================================================================

// Build a vowel frame (PARALLEL topology, RESONATOR type)
inline FilterFrame buildVowelFrame(const VowelFormants& v, double sr) {
    FilterFrame frame;
    frame.topology = Topology::PARALLEL;
    frame.filterType = FilterType::RESONATOR;

    // 4 active formant poles
    frame.poles[0] = ZPlaneDesigner::makePole(v.f1, v.bw, sr, v.w1);
    frame.poles[1] = ZPlaneDesigner::makePole(v.f2, v.bw, sr, v.w2);
    frame.poles[2] = ZPlaneDesigner::makePole(v.f3, v.bw, sr, v.w3);
    frame.poles[3] = ZPlaneDesigner::makePole(v.f4, v.bw, sr, v.w4);

    // Stages 5-7: muted
    frame.poles[4] = ZPlaneDesigner::mutedPole();
    frame.poles[5] = ZPlaneDesigner::mutedPole();
    frame.poles[6] = ZPlaneDesigner::mutedPole();

    // Assign stable IDs
    for (int i = 0; i < NUM_STAGES; ++i) {
        frame.poles[i].id = i;
    }

    return frame;
}

// Build a phaser frame (CASCADE topology, ALLPASS type)
inline FilterFrame buildPhaserFrame(double startFreq, double endFreq,
                                     double bandwidth, double sr,
                                     int stages = 7) {
    FilterFrame frame;
    frame.topology = Topology::CASCADE;
    frame.filterType = FilterType::ALLPASS;

    // Logarithmic frequency spread
    for (int i = 0; i < NUM_STAGES; ++i) {
        if (i < stages) {
            double t = (stages > 1) ? static_cast<double>(i) / (stages - 1) : 0.0;
            double freq = startFreq * std::pow(endFreq / startFreq, t);
            frame.poles[i] = ZPlaneDesigner::makePole(freq, bandwidth, sr);
            frame.poles[i].id = i;
        } else {
            frame.poles[i] = ZPlaneDesigner::mutedPole();
            frame.poles[i].id = i;
        }
    }

    return frame;
}

// Build an acid frame (CASCADE topology, stacked poles at cutoff)
inline FilterFrame buildAcidFrame(double cutoff, double resonance, double sr,
                                   int stages = 4) {
    FilterFrame frame;
    frame.topology = Topology::CASCADE;
    frame.filterType = FilterType::LOWPASS;

    // Convert resonance (0-1) to radius
    double radius = 0.7 + resonance * 0.299;  // 0.7 to 0.999

    for (int i = 0; i < NUM_STAGES; ++i) {
        if (i < stages) {
            Pole p;
            p.freq = cutoff;
            p.radius = radius;
            p.weight = 1.0;
            p.id = i;
            frame.poles[i] = p;
        } else {
            frame.poles[i] = ZPlaneDesigner::mutedPole();
            frame.poles[i].id = i;
        }
    }

    return frame;
}

// Build a custom frame from raw pole specs
inline FilterFrame buildCustomFrame(const Pole* poles, int count,
                                     Topology topology, FilterType type) {
    FilterFrame frame;
    frame.topology = topology;
    frame.filterType = type;

    for (int i = 0; i < NUM_STAGES; ++i) {
        if (i < count) {
            frame.poles[i] = poles[i];
            if (frame.poles[i].id < 0) {
                frame.poles[i].id = i;
            }
        } else {
            frame.poles[i] = ZPlaneDesigner::mutedPole();
            frame.poles[i].id = i;
        }
    }

    return frame;
}

} // namespace field::dsp
