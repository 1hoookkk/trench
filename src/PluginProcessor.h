#pragma once
// =============================================================================
// FIELD — PluginProcessor.h
// =============================================================================
// JUCE AudioProcessor implementation for FIELD Z-Plane filter plugin.
//
// UI Parameters (3-control system):
//   - FIELD (0-1): Position in filter universe, selects vowel pairs
//     - [0.0-0.33): Ah/Ee morph
//     - [0.33-0.66): Oo/Eh morph
//     - [0.66-1.0]: Uh/Ee morph
//   - TENSION (0-1): Maps to drive (0.5-3.0) and postGain (2.0-8.0)
//   - MOTION (0-1): Morph depth between Frame A and Frame B
//   - Test Noise (dev toggle): White noise generator for DSP testing
//
// Architecture:
//   - Dual mono processing (independent L/R ZPlaneCore instances)
//   - Frame selection in updateFramesFromField (cpp:72-100)
//   - Parameter processing in processBlock (cpp:120-175)
//   - Test noise generator using std::mt19937 (cpp:162-163)
//
// Integration:
//   - Uses ZPlaneCore for DSP processing
//   - Uses ZPlaneDesigner for frame generation
//   - APVTS for parameter management
//
// See also:
//   - PluginProcessor.cpp:72-100: FIELD parameter frame mapping
//   - PluginProcessor.cpp:140-148: TENSION parameter mapping
//   - src/dsp/ZPlaneCore.h: DSP engine
//   - src/dsp/ZPlaneDesigner.h: Frame factory
// =============================================================================

#include <juce_audio_processors/juce_audio_processors.h>
#include "dsp/ZPlaneCore.h"
#include "dsp/ZPlaneDesigner.h"
#include <random>

class FIELDProcessor : public juce::AudioProcessor
{
public:
    FIELDProcessor();
    ~FIELDProcessor() override = default;

    void prepareToPlay(double sampleRate, int samplesPerBlock) override;
    void releaseResources() override;
    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;
    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override { return true; }

    const juce::String getName() const override { return JucePlugin_Name; }
    bool acceptsMidi() const override { return false; }
    bool producesMidi() const override { return false; }
    bool isMidiEffect() const override { return false; }
    double getTailLengthSeconds() const override { return 0.5; }

    int getNumPrograms() override { return 1; }
    int getCurrentProgram() override { return 0; }
    void setCurrentProgram(int) override {}
    const juce::String getProgramName(int) override { return {}; }
    void changeProgramName(int, const juce::String&) override {}

    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;

    juce::AudioProcessorValueTreeState& getAPVTS() { return apvts_; }

private:
    juce::AudioProcessorValueTreeState::ParameterLayout createParameterLayout();
    juce::AudioProcessorValueTreeState apvts_;

    // Update frames based on FIELD parameter
    void updateFramesFromField(float fieldValue);

    // DSP
    field::dsp::ZPlaneCore coreL_;
    field::dsp::ZPlaneCore coreR_;
    double sampleRate_ = 44100.0;

    // Frame state
    float lastFieldValue_ = -1.0f;

    // White noise
    std::mt19937 rng_;
    std::uniform_real_distribution<float> noise_{-1.0f, 1.0f};

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FIELDProcessor)
};
