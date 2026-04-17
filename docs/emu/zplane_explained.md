# E-mu Z-plane filter explained

Archived from an E-mu Proteus 2000-era reference manual. This is the
canonical description of the Z-plane filter that TRENCH implements. Do
not paraphrase — this is the ground truth.

## Source text

A Z-plane filter is a filter which can change its function over time.
In a Z-plane filter, we start with two complex filter types and
interpolate between them using a single parameter.

*[Diagram in original manual: a 3D plot with Frequency on the horizontal
axis, Amplitude on the vertical axis, and Morph on the depth axis. Two
filter "frames" labelled `A Filter` and `B Filter` define the two edges
of a 2D surface; the Morph axis interpolates smoothly between them.
Dashed lines mark the A and B endpoints.]*

The Z-plane filter has the unique ability to change its function over time.

Filters A and B represent two different complex filters or "frames."
Changing a single parameter, the Morph, changes many complex filter
parameters simultaneously. Following along the Morph axis you can see
that the filter response smoothly interpolates between the two filters.
This is the essence of the Z-plane filter. Through the use of
interpolation, many complex parameters are condensed into one manageable
entity.

Consider, as an example, the human vocal tract, which is a type of
complex filter or resonator. There are dozens of different muscles
controlling the shape of the vocal tract. When speaking, however, we
don't think of the muscles, we just remember how it feels to form the
vowels. A vowel is really a configuration of many muscles, but we
consider it a single object. In changing from one vowel to another, we
don't need to consider the frequencies of the resonant peaks. You
remember the shape of your mouth for each sound and interpolate between
them.

This Z-plane filter sweep can be controlled by an envelope generator,
an LFO, modulation wheels or pedals, keyboard velocity, key pressure,
and so on. In fact, any of the modulation sources can control the
Z-plane filter.

Because creating the complex filtering is difficult and very time
consuming, we have created 50 different types of filters and installed
them permanently in ROM for your use. You simply select and use the
filters in a manner similar to choosing an instrument. Because there
are so many types of instruments and filters to choose from, the number
of possible permutations is staggering.

## TRENCH implications

- A pill is one Z-plane filter preset: a pair of complex filter frames
  (A and B) plus the Morph axis between them.
- MORPH is driven by any modulation source at playback. TRENCH has no
  custom compiler or time-grid scheduler; the DAW's own mod routing is
  the interface.
- The vocal-tract metaphor is E-mu's own language. TRENCH's phoneme
  pill vocabulary is a faithful translation, not an invented
  abstraction.
- The 50-filter figure in this text is from Proteus 2000-era
  documentation. The Emulator X3 bank has 55 filters; see
  [`emulator_x3_filter_reference.md`](emulator_x3_filter_reference.md).
