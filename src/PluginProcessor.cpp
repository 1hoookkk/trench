// =============================================================================
// FIELD — PluginProcessor.cpp
// =============================================================================
// 3-control Z-Plane engine: FIELD, TENSION, MOTION
// =============================================================================

#include "PluginProcessor.h"
#include "PluginEditor.h"

FIELDProcessor::FIELDProcessor()
    : AudioProcessor(BusesProperties()
          .withInput("Input", juce::AudioChannelSet::stereo(), true)
          .withOutput("Output", juce::AudioChannelSet::stereo(), true))
    , apvts_(*this, nullptr, "Parameters", createParameterLayout())
    , rng_(std::random_device{}())
{
}

juce::AudioProcessorValueTreeState::ParameterLayout FIELDProcessor::createParameterLayout()
{
    std::vector<std::unique_ptr<juce::RangedAudioParameter>> params;

    // FIELD: position in filter universe (0..1)
    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID{"field", 1},
        "Field",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f),
        0.0f
    ));

    // TENSION: Q/drive/compensation intensity (0..1)
    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID{"tension", 1},
        "Tension",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f),
        0.5f
    ));

    // MOTION: morph depth/activity (0..1)
    params.push_back(std::make_unique<juce::AudioParameterFloat>(
        juce::ParameterID{"motion", 1},
        "Motion",
        juce::NormalisableRange<float>(0.0f, 1.0f, 0.001f),
        0.0f
    ));

    // Dev: test noise toggle
    params.push_back(std::make_unique<juce::AudioParameterBool>(
        juce::ParameterID{"testnoise", 1},
        "Test Noise",
        true
    ));

    return {params.begin(), params.end()};
}

void FIELDProcessor::prepareToPlay(double sampleRate, int /*samplesPerBlock*/)
{
    sampleRate_ = sampleRate;

    // Configure both channels
    coreL_.setSampleRate(sampleRate_);
    coreR_.setSampleRate(sampleRate_);

    // Initial frame setup (will be updated by FIELD param)
    updateFramesFromField(0.0f);

    coreL_.reset();
    coreR_.reset();
}

void FIELDProcessor::updateFramesFromField(float fieldValue)
{
    // Simple piecewise mapping: FIELD selects vowel pairs
    // 0.0 = Ah/Ee, 0.5 = Oo/Eh, 1.0 = Uh/Ee
    // Later: replace with procedural generator

    field::dsp::FilterFrame frameA, frameB;

    if (fieldValue < 0.33f) {
        // Cluster A: Ah/Ee (classic vowel morph)
        frameA = field::dsp::buildVowelFrame(field::dsp::VowelFormants::Ah(), sampleRate_);
        frameB = field::dsp::buildVowelFrame(field::dsp::VowelFormants::Ee(), sampleRate_);
    } else if (fieldValue < 0.66f) {
        // Cluster B: Oo/Eh
        frameA = field::dsp::buildVowelFrame(field::dsp::VowelFormants::Oo(), sampleRate_);
        frameB = field::dsp::buildVowelFrame(field::dsp::VowelFormants::Eh(), sampleRate_);
    } else {
        // Cluster C: Uh/Ee
        frameA = field::dsp::buildVowelFrame(field::dsp::VowelFormants::Uh(), sampleRate_);
        frameB = field::dsp::buildVowelFrame(field::dsp::VowelFormants::Ee(), sampleRate_);
    }

    coreL_.setFrameA(frameA);
    coreL_.setFrameB(frameB);
    coreR_.setFrameA(frameA);
    coreR_.setFrameB(frameB);

    lastFieldValue_ = fieldValue;
}

void FIELDProcessor::releaseResources()
{
    coreL_.reset();
    coreR_.reset();
}

bool FIELDProcessor::isBusesLayoutSupported(const BusesLayout& layouts) const
{
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono()
        && layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo())
        return false;

    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet())
        return false;

    return true;
}

void FIELDProcessor::processBlock(juce::AudioBuffer<float>& buffer,
                                   juce::MidiBuffer& /*midiMessages*/)
{
    juce::ScopedNoDenormals noDenormals;

    const int numSamples = buffer.getNumSamples();
    const int numChannels = buffer.getNumChannels();

    // Read params
    float fieldValue = apvts_.getRawParameterValue("field")->load();
    float tension = apvts_.getRawParameterValue("tension")->load();
    float motion = apvts_.getRawParameterValue("motion")->load();
    bool useNoise = apvts_.getRawParameterValue("testnoise")->load() > 0.5f;

    // Update frames if FIELD changed significantly
    if (std::abs(fieldValue - lastFieldValue_) > 0.01f) {
        updateFramesFromField(fieldValue);
    }

    // TENSION maps to drive + postGain
    // Low tension = gentle, high tension = scream
    float drive = 0.5f + tension * 2.5f;       // 0.5 to 3.0
    float postGain = 2.0f + tension * 6.0f;    // 2.0 to 8.0

    coreL_.setDrive(drive);
    coreL_.setPostGain(postGain);
    coreR_.setDrive(drive);
    coreR_.setPostGain(postGain);

    // MOTION controls morph depth
    // At motion=0, morph stays at 0 (Frame A only)
    // At motion=1, morph sweeps full range
    float morphValue = motion;  // For now: direct mapping
    // Later: add LFO/envelope modulation

    coreL_.setMorph(morphValue);
    coreR_.setMorph(morphValue);

    // Process left channel
    if (numChannels >= 1) {
        float* L = buffer.getWritePointer(0);
        for (int i = 0; i < numSamples; ++i) {
            float input = useNoise ? noise_(rng_) : L[i];
            L[i] = coreL_.processSample(input);
        }
    }

    // Process right channel
    if (numChannels >= 2) {
        float* R = buffer.getWritePointer(1);
        for (int i = 0; i < numSamples; ++i) {
            float input = useNoise ? noise_(rng_) : R[i];
            R[i] = coreR_.processSample(input);
        }
    }
}

juce::AudioProcessorEditor* FIELDProcessor::createEditor()
{
    return new FIELDEditor(*this);
}

void FIELDProcessor::getStateInformation(juce::MemoryBlock& destData)
{
    auto state = apvts_.copyState();
    std::unique_ptr<juce::XmlElement> xml(state.createXml());
    copyXmlToBinary(*xml, destData);
}

void FIELDProcessor::setStateInformation(const void* data, int sizeInBytes)
{
    std::unique_ptr<juce::XmlElement> xml(getXmlFromBinary(data, sizeInBytes));
    if (xml && xml->hasTagName(apvts_.state.getType()))
        apvts_.replaceState(juce::ValueTree::fromXml(*xml));
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter()
{
    return new FIELDProcessor();
}
