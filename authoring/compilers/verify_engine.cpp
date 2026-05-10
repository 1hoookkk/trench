// =============================================================================
// verify_engine.cpp — Offline verification harness for FIELD DSP
// =============================================================================
// Renders test signals, exports WAV + stats for comparison against references
// Build: cl /EHsc /O2 /I../src verify_engine.cpp /Fe:verify_engine.exe
// =============================================================================

#include <cmath>
#include <cstdio>
#include <cstdint>
#include <vector>
#include <random>
#include <algorithm>

// Include the DSP headers directly (header-only)
#include "../src/dsp/ZPlaneCore.h"
#include "../src/dsp/ZPlaneDesigner.h"

using namespace field::dsp;

// =============================================================================
// WAV Writer (minimal, 32-bit float)
// =============================================================================

void writeWavFloat(const char* filename, const float* data, size_t numSamples, double sampleRate) {
    FILE* f = fopen(filename, "wb");
    if (!f) { printf("Failed to open %s\n", filename); return; }

    uint32_t dataSize = static_cast<uint32_t>(numSamples * sizeof(float));
    uint32_t fileSize = 36 + dataSize;

    // RIFF header
    fwrite("RIFF", 1, 4, f);
    fwrite(&fileSize, 4, 1, f);
    fwrite("WAVE", 1, 4, f);

    // fmt chunk (IEEE float)
    fwrite("fmt ", 1, 4, f);
    uint32_t fmtSize = 16;
    uint16_t audioFormat = 3;  // IEEE float
    uint16_t numChannels = 1;
    uint32_t sr = static_cast<uint32_t>(sampleRate);
    uint32_t byteRate = sr * sizeof(float);
    uint16_t blockAlign = sizeof(float);
    uint16_t bitsPerSample = 32;

    fwrite(&fmtSize, 4, 1, f);
    fwrite(&audioFormat, 2, 1, f);
    fwrite(&numChannels, 2, 1, f);
    fwrite(&sr, 4, 1, f);
    fwrite(&byteRate, 4, 1, f);
    fwrite(&blockAlign, 2, 1, f);
    fwrite(&bitsPerSample, 2, 1, f);

    // data chunk
    fwrite("data", 1, 4, f);
    fwrite(&dataSize, 4, 1, f);
    fwrite(data, sizeof(float), numSamples, f);

    fclose(f);
    printf("Wrote %s (%zu samples)\n", filename, numSamples);
}

// =============================================================================
// Signal Generators
// =============================================================================

void generateWhiteNoise(float* buffer, size_t numSamples, uint32_t seed = 12345) {
    std::mt19937 rng(seed);
    std::uniform_real_distribution<float> dist(-1.0f, 1.0f);
    for (size_t i = 0; i < numSamples; ++i) {
        buffer[i] = dist(rng);
    }
}

void generateImpulse(float* buffer, size_t numSamples) {
    std::fill(buffer, buffer + numSamples, 0.0f);
    if (numSamples > 0) buffer[0] = 1.0f;
}

void generateSaw(float* buffer, size_t numSamples, double freq, double sampleRate) {
    double phase = 0.0;
    double inc = freq / sampleRate;
    for (size_t i = 0; i < numSamples; ++i) {
        buffer[i] = static_cast<float>(2.0 * phase - 1.0);
        phase += inc;
        if (phase >= 1.0) phase -= 1.0;
    }
}

// =============================================================================
// Statistics
// =============================================================================

struct AudioStats {
    float peak;
    float rms;
    float crest;  // peak/rms in dB
};

AudioStats computeStats(const float* buffer, size_t numSamples) {
    AudioStats s{};
    double sumSq = 0.0;
    for (size_t i = 0; i < numSamples; ++i) {
        float absVal = std::abs(buffer[i]);
        if (absVal > s.peak) s.peak = absVal;
        sumSq += buffer[i] * buffer[i];
    }
    s.rms = static_cast<float>(std::sqrt(sumSq / numSamples));
    s.crest = (s.rms > 1e-12f) ? 20.0f * std::log10(s.peak / s.rms) : 0.0f;
    return s;
}

void printStats(const char* name, const AudioStats& s) {
    printf("  %-20s  Peak: %.4f (%.1f dB)  RMS: %.4f (%.1f dB)  Crest: %.1f dB\n",
           name, s.peak, 20.0 * std::log10(s.peak + 1e-12),
           s.rms, 20.0 * std::log10(s.rms + 1e-12), s.crest);
}

// =============================================================================
// Frequency Response (via impulse)
// =============================================================================

void writeFrequencyResponse(const char* filename, const float* impulse, size_t numSamples, double sampleRate) {
    FILE* f = fopen(filename, "w");
    if (!f) return;

    fprintf(f, "freq_hz,magnitude_db\n");

    // Simple DFT for key frequencies (not FFT, but sufficient for verification)
    const double freqs[] = {50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000,
                            1200, 1500, 2000, 2500, 3000, 4000, 5000, 6000, 8000,
                            10000, 12000, 15000, 18000, 20000};

    for (double freq : freqs) {
        if (freq >= sampleRate / 2) break;

        double realSum = 0.0, imagSum = 0.0;
        double omega = 2.0 * PI * freq / sampleRate;

        for (size_t n = 0; n < numSamples; ++n) {
            realSum += impulse[n] * std::cos(omega * n);
            imagSum -= impulse[n] * std::sin(omega * n);
        }

        double mag = std::sqrt(realSum * realSum + imagSum * imagSum);
        double db = 20.0 * std::log10(mag + 1e-12);

        fprintf(f, "%.1f,%.2f\n", freq, db);
    }

    fclose(f);
    printf("Wrote %s\n", filename);
}

// =============================================================================
// Test Cases
// =============================================================================

void testVowelAh(double sampleRate) {
    printf("\n=== Test: Vowel Ah (white noise input) ===\n");

    ZPlaneCore core;
    core.setSampleRate(sampleRate);

    auto frameAh = buildVowelFrame(VowelFormants::Ah(), sampleRate);
    core.setFrameA(frameAh);
    core.setFrameB(frameAh);  // No morph
    core.setMorph(0.0f);
    core.reset();

    const size_t numSamples = static_cast<size_t>(sampleRate * 5.0);  // 5 seconds
    std::vector<float> input(numSamples);
    std::vector<float> output(numSamples);

    generateWhiteNoise(input.data(), numSamples);

    for (size_t i = 0; i < numSamples; ++i) {
        output[i] = core.processSample(input[i]);
    }

    printStats("Input (noise)", computeStats(input.data(), numSamples));
    printStats("Output (Ah)", computeStats(output.data(), numSamples));

    writeWavFloat("verify_vowel_ah_noise_in.wav", input.data(), numSamples, sampleRate);
    writeWavFloat("verify_vowel_ah_out.wav", output.data(), numSamples, sampleRate);

    // Impulse response
    std::vector<float> impulseIn(8192, 0.0f);
    std::vector<float> impulseOut(8192);
    impulseIn[0] = 1.0f;
    core.reset();
    for (size_t i = 0; i < 8192; ++i) {
        impulseOut[i] = core.processSample(impulseIn[i]);
    }
    writeFrequencyResponse("verify_vowel_ah_fr.csv", impulseOut.data(), 8192, sampleRate);
}

void testMorphAhToEe(double sampleRate) {
    printf("\n=== Test: Morph Ah -> Ee (saw input) ===\n");

    ZPlaneCore core;
    core.setSampleRate(sampleRate);

    auto frameAh = buildVowelFrame(VowelFormants::Ah(), sampleRate);
    auto frameEe = buildVowelFrame(VowelFormants::Ee(), sampleRate);
    core.setFrameA(frameAh);
    core.setFrameB(frameEe);
    core.reset();

    const size_t numSamples = static_cast<size_t>(sampleRate * 10.0);  // 10 seconds
    std::vector<float> input(numSamples);
    std::vector<float> output(numSamples);

    generateSaw(input.data(), numSamples, 110.0, sampleRate);  // A2

    // Morph from 0 to 1 over duration
    for (size_t i = 0; i < numSamples; ++i) {
        float t = static_cast<float>(i) / static_cast<float>(numSamples - 1);
        core.setMorph(t);
        output[i] = core.processSample(input[i]);
    }

    printStats("Input (saw 110Hz)", computeStats(input.data(), numSamples));
    printStats("Output (Ah->Ee)", computeStats(output.data(), numSamples));

    writeWavFloat("verify_morph_ah_ee_saw_in.wav", input.data(), numSamples, sampleRate);
    writeWavFloat("verify_morph_ah_ee_out.wav", output.data(), numSamples, sampleRate);
}

void testPhaser(double sampleRate) {
    printf("\n=== Test: Phaser (saw input) ===\n");

    ZPlaneCore core;
    core.setSampleRate(sampleRate);

    auto frameA = buildPhaserFrame(200.0, 4000.0, 100.0, sampleRate, 6);
    auto frameB = buildPhaserFrame(400.0, 8000.0, 100.0, sampleRate, 6);
    core.setFrameA(frameA);
    core.setFrameB(frameB);
    core.setPostGain(1.0f);  // Phasers don't need huge makeup
    core.reset();

    const size_t numSamples = static_cast<size_t>(sampleRate * 5.0);
    std::vector<float> input(numSamples);
    std::vector<float> output(numSamples);

    generateSaw(input.data(), numSamples, 220.0, sampleRate);

    // Sweep morph
    for (size_t i = 0; i < numSamples; ++i) {
        float t = 0.5f + 0.5f * std::sin(2.0 * PI * 0.2 * i / sampleRate);
        core.setMorph(t);
        output[i] = core.processSample(input[i]);
    }

    printStats("Input (saw 220Hz)", computeStats(input.data(), numSamples));
    printStats("Output (phaser)", computeStats(output.data(), numSamples));

    writeWavFloat("verify_phaser_out.wav", output.data(), numSamples, sampleRate);
}

void printCoefficients(double sampleRate) {
    printf("\n=== Coefficient Dump (for ROM comparison) ===\n");

    auto frameAh = buildVowelFrame(VowelFormants::Ah(), sampleRate);

    printf("Vowel Ah @ %.0f Hz:\n", sampleRate);
    printf("  Stage  Freq(Hz)    Radius      a1            a2            b0            b2\n");

    for (int i = 0; i < NUM_STAGES; ++i) {
        const Pole& p = frameAh.poles[i];
        if (!p.isActive()) {
            printf("  [%d]    (muted)\n", i);
            continue;
        }

        double theta = 2.0 * PI * p.freq / sampleRate;
        double R = p.radius;
        double a1 = -2.0 * R * std::cos(theta);
        double a2 = R * R;
        double gain = 1.0 - R;
        double b0 = gain;
        double b2 = -gain;

        printf("  [%d]    %7.1f    %.6f    %+.8f    %+.8f    %+.8f    %+.8f\n",
               i, p.freq, R, a1, a2, b0, b2);
    }
}

// =============================================================================
// Main
// =============================================================================

int main() {
    printf("FIELD DSP Verification Harness\n");
    printf("==============================\n");

    const double sampleRate = 48000.0;

    printCoefficients(sampleRate);
    testVowelAh(sampleRate);
    testMorphAhToEe(sampleRate);
    testPhaser(sampleRate);

    printf("\n=== Done ===\n");
    printf("Compare WAV outputs against X3 reference captures.\n");
    printf("Compare FR CSV against analyzer traces.\n");

    return 0;
}
