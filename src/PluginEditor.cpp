// =============================================================================
// FIELD — PluginEditor.cpp
// =============================================================================
// 3-control UI: FIELD, TENSION, MOTION
// =============================================================================

#include "PluginEditor.h"

FIELDEditor::FIELDEditor(FIELDProcessor& processor)
    : AudioProcessorEditor(&processor)
    , processor_(processor)
{
    setSize(320, 220);

    auto setupSlider = [this](juce::Slider& slider, juce::Label& label, const juce::String& name) {
        slider.setSliderStyle(juce::Slider::LinearHorizontal);
        slider.setTextBoxStyle(juce::Slider::TextBoxRight, false, 45, 20);
        slider.setColour(juce::Slider::thumbColourId, juce::Colour(0xFF70B0A0));
        slider.setColour(juce::Slider::trackColourId, juce::Colour(0xFF404040));
        slider.setColour(juce::Slider::backgroundColourId, juce::Colour(0xFF202020));
        addAndMakeVisible(slider);

        label.setText(name, juce::dontSendNotification);
        label.setJustificationType(juce::Justification::left);
        label.setColour(juce::Label::textColourId, juce::Colour(0xFF808080));
        addAndMakeVisible(label);
    };

    setupSlider(fieldSlider_, fieldLabel_, "FIELD");
    setupSlider(tensionSlider_, tensionLabel_, "TENSION");
    setupSlider(motionSlider_, motionLabel_, "MOTION");

    // Test noise toggle
    noiseToggle_.setButtonText("Test Noise");
    noiseToggle_.setColour(juce::ToggleButton::textColourId, juce::Colour(0xFF606060));
    noiseToggle_.setColour(juce::ToggleButton::tickColourId, juce::Colour(0xFF70B0A0));
    addAndMakeVisible(noiseToggle_);

    // Attach to parameters
    auto& apvts = processor_.getAPVTS();
    fieldAttach_ = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>(
        apvts, "field", fieldSlider_);
    tensionAttach_ = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>(
        apvts, "tension", tensionSlider_);
    motionAttach_ = std::make_unique<juce::AudioProcessorValueTreeState::SliderAttachment>(
        apvts, "motion", motionSlider_);
    noiseAttach_ = std::make_unique<juce::AudioProcessorValueTreeState::ButtonAttachment>(
        apvts, "testnoise", noiseToggle_);
}

void FIELDEditor::paint(juce::Graphics& g)
{
    g.fillAll(juce::Colour(0xFF111111));

    // Title
    g.setColour(juce::Colour(0xFF70B0A0));
    g.setFont(juce::FontOptions(juce::Font::getDefaultMonospacedFontName(), 20.0f, juce::Font::bold));
    g.drawText("FIELD", 10, 10, getWidth() - 20, 28, juce::Justification::centred);

    // Subtitle
    g.setColour(juce::Colour(0xFF505050));
    g.setFont(juce::FontOptions(juce::Font::getDefaultMonospacedFontName(), 10.0f, juce::Font::plain));
    g.drawText("Z-Plane Filter Generator", 10, 38, getWidth() - 20, 14, juce::Justification::centred);

    // Hairline separator
    g.setColour(juce::Colour(0xFF303030));
    g.drawHorizontalLine(55, 15, getWidth() - 15);
}

void FIELDEditor::resized()
{
    const int margin = 15;
    const int labelW = 60;
    const int sliderH = 28;
    const int gap = 8;
    int y = 65;

    auto layoutRow = [&](juce::Label& label, juce::Slider& slider) {
        label.setBounds(margin, y, labelW, sliderH);
        slider.setBounds(margin + labelW, y, getWidth() - margin * 2 - labelW, sliderH);
        y += sliderH + gap;
    };

    layoutRow(fieldLabel_, fieldSlider_);
    layoutRow(tensionLabel_, tensionSlider_);
    layoutRow(motionLabel_, motionSlider_);

    y += 5;
    noiseToggle_.setBounds(margin, y, getWidth() - margin * 2, 25);
}
