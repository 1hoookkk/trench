#pragma once
// =============================================================================
// FIELD — PluginEditor.h
// =============================================================================
// 3-control UI: FIELD, TENSION, MOTION
// =============================================================================

#include <juce_audio_processors/juce_audio_processors.h>
#include "PluginProcessor.h"

class FIELDEditor : public juce::AudioProcessorEditor
{
public:
    explicit FIELDEditor(FIELDProcessor& processor);
    ~FIELDEditor() override = default;

    void paint(juce::Graphics&) override;
    void resized() override;

private:
    FIELDProcessor& processor_;

    // 3 main controls
    juce::Slider fieldSlider_;
    juce::Slider tensionSlider_;
    juce::Slider motionSlider_;

    juce::Label fieldLabel_;
    juce::Label tensionLabel_;
    juce::Label motionLabel_;

    // Dev toggle
    juce::ToggleButton noiseToggle_;

    // Attachments
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> fieldAttach_;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> tensionAttach_;
    std::unique_ptr<juce::AudioProcessorValueTreeState::SliderAttachment> motionAttach_;
    std::unique_ptr<juce::AudioProcessorValueTreeState::ButtonAttachment> noiseAttach_;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(FIELDEditor)
};
