# fourier-synthesis by Joseph Percival
This repository implements an experimental module for use in VCVRack.

Fourier Synthesis asks the question: what if the Fourier transform used sawtooth waves instead of sinusoids? This module takes a stereo input and passes it through a discrete Fourier transform (DFT), resulting in a real-time frequency analysis of the input signal (presented on the display). Next, this frequency data is modified according to the selected waveform type, so that when the data is passed back through an inverse DFT, the output audio is reconstructed from a series of that particular waveform (e.g. a series of sawtooths) rather than a series of sinusoids as per the regular algorithm. The result is that even a very simple audio input such as a square wave can very quickly become distorted as each harmonic in the square wave generates its own harmonic. 

The module makes use of the FFTW library, which can handle any integer buffer size selected by the buffer size parameter. Additional parameters include sample rate reduction, waveform shape, and number of harmonics, giving plenty of options to experiment with.

I recommend starting off with a small buffer size and slowly increasing it as you experiment with the other parameters. There is also the option of daisy-chaining this module; you can also find some interesting effects with just one module by plugging the left output into the right input or vice versa (connecting left and left or right and right will cause a feedback loop).

## Setup Instructions

1. Clone this repository into the relevant plugins folder according to the instructions listed [here](https://vcvrack.com/manual/Building).
1. Run `./setup.sh` to download and compile the FFTW library.
2. Build the module with `make`.
