# fourier-synthesis by Joseph Percival
This repository implements a module for use in VCVRack. 

The fourier algorithm is used to model a waveform which is input to the module. This model is an approximation of the original input, and is made up of a series of sine waves of varying frequency, phase and amplitude. 

This modelling can be done in real time by sampling the input waveform thousands of times per second. Each tiny sample can be recorded into an input buffer, which the fourier algorithm can then be applied to to generate a series of sine waves which approximate the sample. This series of sine waves is then combined into a single waveform which is then output by the module.

The motivation for this is that by using the fourier transform we are essentially performing real-time additive synthesis to model a waveform, meaning each of the individual component sine waves of the input waveform can be modified without affecting the rest of the sound. This module can also be treated as an effects module - the fourier transform is only a perfect model with an infinite number of sine waves, so the waveform output will have some distortion or audio artifacts when compared to the original. The exact results of this will be apparent once the module is functional.
