# ==============================================================================
# FIELD — Z-Plane Generator Plugin
# ==============================================================================
# "We don't emulate the ROM. We emulate the engineer."
# ==============================================================================

cmake_minimum_required(VERSION 3.22)

# Version
set(PROJECT_VERSION 0.0.1)

# Plugin identifiers
set(PLUGIN_NAME "FIELD")
set(PLUGIN_CODE "Fld1")
set(MANUFACTURER_NAME "Trench")
set(MANUFACTURER_CODE "Trch")

project(${PLUGIN_NAME} VERSION ${PROJECT_VERSION})

# C++20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ==============================================================================
# JUCE via FetchContent (Pamplejuce pattern)
# ==============================================================================

include(FetchContent)

FetchContent_Declare(
    JUCE
    GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
    GIT_TAG 8.0.10
    GIT_SHALLOW TRUE
    GIT_PROGRESS TRUE
)
FetchContent_MakeAvailable(JUCE)

# Build options
option(FIELD_COPY_PLUGIN_AFTER_BUILD "Copy plugin to system folder after build" ON)

# Platform-specific formats
if(APPLE)
    set(FIELD_FORMATS VST3 AU Standalone)
else()
    set(FIELD_FORMATS VST3 Standalone)
endif()

# ==============================================================================
# Plugin Target
# ==============================================================================

juce_add_plugin(${PLUGIN_NAME}
    VERSION ${PROJECT_VERSION}
    COMPANY_NAME "${MANUFACTURER_NAME}"
    BUNDLE_ID "com.trench.field"
    IS_SYNTH FALSE
    NEEDS_MIDI_INPUT FALSE
    NEEDS_MIDI_OUTPUT FALSE
    IS_MIDI_EFFECT FALSE
    EDITOR_WANTS_KEYBOARD_FOCUS FALSE
    COPY_PLUGIN_AFTER_BUILD ${FIELD_COPY_PLUGIN_AFTER_BUILD}
    PLUGIN_MANUFACTURER_CODE ${MANUFACTURER_CODE}
    PLUGIN_CODE ${PLUGIN_CODE}
    FORMATS ${FIELD_FORMATS}
    PRODUCT_NAME "${PLUGIN_NAME}"
)

# ==============================================================================
# Sources
# ==============================================================================

target_sources(${PLUGIN_NAME}
    PRIVATE
        src/PluginProcessor.cpp
        src/PluginEditor.cpp
)

# ==============================================================================
# Dependencies
# ==============================================================================

target_compile_definitions(${PLUGIN_NAME}
    PUBLIC
        JUCE_WEB_BROWSER=0
        JUCE_USE_CURL=0
        JUCE_VST3_CAN_REPLACE_VST2=0
)

target_include_directories(${PLUGIN_NAME}
    PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/src
        ${CMAKE_CURRENT_SOURCE_DIR}/src/dsp
)

target_link_libraries(${PLUGIN_NAME}
    PRIVATE
        juce::juce_audio_utils
        juce::juce_dsp
    PUBLIC
        juce::juce_recommended_config_flags
        juce::juce_recommended_lto_flags
        juce::juce_recommended_warning_flags
)
