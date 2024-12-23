#include "fourier-synthesis-plugin.hpp"
#include <fftw3.h>
#include <vector>
#include <iostream>

struct FourierSynthesis : Module {
    int bufferSize;
    int sampleRate;
    float waveformType;
    int numHarmonics;
    fftw_plan forward_plan;
    fftw_plan backward_plan;
    double* real_in;
    fftw_complex* freq_out;
    double* real_out;
    int bufferIndex;
    int sampleRateIndex;
    std::vector<double> freqMagnitudes;

    enum ParamIds {
        BUFFER_PARAM,
        SAMPLE_RATE_PARAM,
        WAVEFORM_PARAM,
        HARMONICS_PARAM,
        NUM_PARAMS
    };

    enum InputIds {
        INPUT_LEFT,
        INPUT_RIGHT,
        SAMPLE_RATE_MOD,
        WAVEFORM_MOD,
        HARMONICS_MOD,
        NUM_INPUTS
    };

    enum OutputIds {
        OUTPUT_LEFT,
        OUTPUT_RIGHT,
        NUM_OUTPUTS
    };
    enum LightId {
        NUM_LIGHTS
    };

    FourierSynthesis() {
        config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
        configParam(BUFFER_PARAM,1.f,10000.f,1.f,"Buffer Size");
		getParamQuantity(BUFFER_PARAM)->snapEnabled = true;
        configParam(SAMPLE_RATE_PARAM,0.f,400.f,1.f,"Sample rate reduction");
		getParamQuantity(SAMPLE_RATE_PARAM)->snapEnabled = true;
        configParam(WAVEFORM_PARAM,0.f,2.f,1.f,"Waveform Type");
        configParam(HARMONICS_PARAM,1.f,100.f,10.f,"Number of Harmonics");
		getParamQuantity(HARMONICS_PARAM)->snapEnabled = true;
        //initialise params & buffers
        bufferSize = params[BUFFER_PARAM].getValue();
        sampleRate = params[SAMPLE_RATE_PARAM].getValue();
        waveformType = params[WAVEFORM_PARAM].getValue();
        numHarmonics = params[HARMONICS_PARAM].getValue();
        bufferIndex = 0;
        sampleRateIndex = 0;
        initialiseResources();
    }

        ~FourierSynthesis() {
        if (real_in) fftw_free(real_in);
        if (real_out) fftw_free(real_out);
        if (freq_out) fftw_free(freq_out);
        if (forward_plan) fftw_destroy_plan(forward_plan);
        if (backward_plan) fftw_destroy_plan(backward_plan);
    }

    void initialiseResources() {
        // free resources
        if (real_in) fftw_free(real_in);
        if (freq_out) fftw_free(freq_out);
        if (real_out) fftw_free(real_out);
        if (forward_plan) fftw_destroy_plan(forward_plan);
        if (backward_plan) fftw_destroy_plan(backward_plan);

        // allocate buffers
        real_in = fftw_alloc_real(bufferSize);
        freq_out = fftw_alloc_complex(bufferSize / 2 + 1);
        real_out = fftw_alloc_real(bufferSize);

        // set real_out array to zero to prevent the module
        // from outputting undefined array content
        memset(real_out, 0, bufferSize * sizeof(double));

        // generate plans
        forward_plan = fftw_plan_dft_r2c_1d(bufferSize, real_in, freq_out, FFTW_ESTIMATE);
        backward_plan = fftw_plan_dft_c2r_1d(bufferSize, freq_out, real_out, FFTW_ESTIMATE);
    }

    void process(const ProcessArgs& args) override {

        if (paramsModified()){
            // reevaluate regular params
            sampleRate = params[SAMPLE_RATE_PARAM].getValue();
            waveformType = params[WAVEFORM_PARAM].getValue();
            numHarmonics = params[HARMONICS_PARAM].getValue();
        }

        if (fftParamsModified()){
            // reevaluate fft params, recalculating plans
            bufferSize = params[BUFFER_PARAM].getValue();
            bufferIndex = 0;
            initialiseResources();
        }

        if (sampleRateIndex < sampleRate) {
            sampleRateIndex++;
        } else {
            sampleRateIndex = 0;
            if (bufferIndex < bufferSize) {
                // stream data from the input to the buffer
                real_in[bufferIndex] = static_cast<double>(inputs[INPUT_LEFT].getVoltage());
                // simultaneously output the data from the previous buffer
                // divide by bufferSize since the output of FFTW is scaled by the size of the input buffer

                outputs[OUTPUT_LEFT].setVoltage(real_out[bufferIndex] / bufferSize);
                outputs[OUTPUT_RIGHT].setVoltage(inputs[INPUT_RIGHT].getVoltage());
                bufferIndex++;
            } else {
                // make sure you don't miss a sample here
                bufferIndex = 0;
                fftw_execute(forward_plan);

                // modify frequency domain with custom waveforms based on chosen wavetable
                applyCustomWaveform(waveformType, bufferSize, freq_out);

                // compute magnitudes for the frequency domain (display)
                freqMagnitudes.resize(bufferSize / 2 + 1);
                for (int i = 0; i < bufferSize / 2 + 1; i++) {
                    freqMagnitudes[i] = sqrt(freq_out[i][0] * freq_out[i][0] + freq_out[i][1] * freq_out[i][1]);
                }

                fftw_execute(backward_plan);
            }
        }
    }

    bool paramsModified() {
        return params[SAMPLE_RATE_PARAM].getValue() != sampleRate ||
               params[WAVEFORM_PARAM].getValue() != waveformType ||
               params[HARMONICS_PARAM].getValue() != numHarmonics;
    }

    bool fftParamsModified() {
        return params[BUFFER_PARAM].getValue() != bufferSize;
    }

    void applyCustomWaveform(double waveformType, int bufferSize, fftw_complex* freq_out) {
        fftw_complex* temp_freq_out = (fftw_complex*) fftw_malloc((bufferSize / 2 + 1) * sizeof(fftw_complex));
        memset(temp_freq_out, 0, (bufferSize / 2 + 1) * sizeof(fftw_complex));
        // uncomment for testing
        // memset(freq_out, 0, (bufferSize / 2 + 1) * sizeof(fftw_complex));
        // freq_out[1][0] = 100;
        // freq_out[1][1] = 50;
        // freq_out[10][0] = 100;
        // freq_out[10][1] = 0;
        // copy DC component
        temp_freq_out[0][0] = freq_out[0][0]; 
        temp_freq_out[0][1] = freq_out[0][1];
        for (int bin = 1; bin < bufferSize / 2 + 1; ++bin) {
            // extract the magnitude and phase of the current bin
            double magnitude = sqrt(freq_out[bin][0] * freq_out[bin][0] + freq_out[bin][1] * freq_out[bin][1]);
            double phase = atan2(freq_out[bin][1], freq_out[bin][0]);
            
            // generate harmonics based on the waveform type
            for (int harmonic = 1; (bin * harmonic < bufferSize / 2 + 1) && (harmonic <= numHarmonics); ++harmonic) {
                int targetBin = bin * harmonic;
                double harmonicMagnitude = magnitude;
                float harmonicCoefficient = waveformType;
                double r;
                double i;

                if (waveformType == 0) { // sine (no harmonics)
                    if (harmonic != 1) continue;
                    r = harmonicMagnitude * cos(phase * harmonic);
                    i = harmonicMagnitude * sin(phase * harmonic);
                } else {
                    if (waveformType <= 1) { // sawtooth
                        harmonicMagnitude /= harmonic;
                    } else if (waveformType <= 2) { // square
                        harmonicMagnitude /= harmonic;
                        if (harmonic % 2 == 1) {
                            harmonicCoefficient = 1;
                        } else {
                            harmonicCoefficient -= 1; 
                            harmonicCoefficient = 1 - harmonicCoefficient;
                        }

                    }
                    if (harmonic == 1) {
                        harmonicCoefficient = 1;
                    }
                    // r = harmonicMagnitude * sin(phase * harmonic);
                    // i = harmonicMagnitude * -cos(phase * harmonic);
                    r = harmonicCoefficient * harmonicMagnitude * cos(phase * harmonic-M_PI_2); // theoretically this 90° phase shift is not needed (-π)
                    i = harmonicCoefficient * harmonicMagnitude * sin(phase * harmonic-M_PI_2); // however when removed the harmonics become misaligned
                    // the next steps would be to compare the output with respect to the DC component
                    // i.e. after applyCustomWaveform, take the fftw_execute(backward_plan) output
                    // and feed it into the forward plan again, comparing both sets of frequency bins -
                    // the phases & magnitudes of the second set are ideally the same as the first
                    // but since we are manually updating the bins in applyCustomWaveform and the DC component
                    // is not updated, the inverse transform may not be accurate.
                }

                // accumulate into the target bin
                temp_freq_out[targetBin][0] += r; // real part
                temp_freq_out[targetBin][1] += i; // imaginary part
            }
        }
        // copy the modified bins back into freq_out
        memcpy(freq_out, temp_freq_out, (bufferSize / 2 + 1) * sizeof(fftw_complex));
        fftw_free(temp_freq_out);
    }

};

struct FrequencyDisplay : Widget {
    FourierSynthesis* module;
    std::vector<double>* freqData = nullptr; // pointer to frequency data
    int numBins = 0;                         // number of bins in display

    void setNumBins(int bins) {
        numBins = bins;
    }

    // map values between 0 and 1
    void scale(std::vector<double>& data) {
        double min_value = *std::min_element(data.begin(), data.end());
        double max_value = *std::max_element(data.begin(), data.end());

        // handle division by 0 error
        if (max_value - min_value == 0) return;

        // normalize each value in the array
        for (auto& value : data) {
            value = (value - min_value) / (max_value - min_value);
        }
    }

    void draw(const DrawArgs& args) override {
        // ensure frequency data is available
        if (!freqData || freqData->empty()) return;

        // get widget size
        NVGcontext* vg = args.vg;
        float width = box.size.x;
        float height = box.size.y;

        // update number of bins
        numBins = (*freqData).size();
        (*freqData).resize(numBins, 0.0f);

        // draw background (temporary for positioning)
        // nvgBeginPath(vg);
        // nvgRect(vg, 0, 0, width, height);
        // nvgFillColor(vg, nvgRGB(20, 20, 20));
        // nvgFill(vg);

        // normalise data
        scale(*freqData);

        // draw frequency spectrum
        nvgBeginPath(vg);
        for (int i = 0; i < numBins; i++) {
            float x = (float)i / numBins * width;
            float barHeight = (*freqData)[i] * height; // scale frequency magnitude
            nvgRect(vg, x, height - barHeight, width / numBins, barHeight);
        }
        // nvgFillColor(vg, nvgRGB(0, 200, 255)); //cyan
        // nvgFillColor(vg, nvgRGB(0, 0, 0)); //black
        nvgFillColor(vg, nvgRGB(230, 233, 169)); //LED white
        // nvgFillColor(vg, nvgRGB(255, 255, 255)); //white
        nvgFill(vg);
    }
};


struct FourierSynthesisWidget : ModuleWidget {
    FrequencyDisplay* display;

    FourierSynthesisWidget(FourierSynthesis* module) {
        setModule(module);
        setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/fourier_bg.svg")));

        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
        // addChild(createWidget<ScrewSilver>(Vec(0, 0)));
        // addChild(createWidget<ScrewSilver>(Vec(box.size.x - RACK_GRID_WIDTH, 0)));
        // addChild(createWidget<ScrewSilver>(Vec(0, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
        // addChild(createWidget<ScrewSilver>(Vec(box.size.x - RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

        addInput(createInput<PJ301MPort>(Vec(18,329), module, FourierSynthesis::INPUT_LEFT));
        addInput(createInput<PJ301MPort>(Vec(47,329), module, FourierSynthesis::INPUT_RIGHT));
        addInput(createInput<PJ301MPort>(Vec(18,10), module, FourierSynthesis::SAMPLE_RATE_MOD));
        addInput(createInput<PJ301MPort>(Vec(47,20), module, FourierSynthesis::WAVEFORM_MOD));
        addInput(createInput<PJ301MPort>(Vec(47,30), module, FourierSynthesis::HARMONICS_MOD));

        addParam(createParam<RoundLargeBlackKnob>(Vec(34,197), module, FourierSynthesis::BUFFER_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(57,235), module, FourierSynthesis::SAMPLE_RATE_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(135,197), module, FourierSynthesis::WAVEFORM_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(158,235), module, FourierSynthesis::HARMONICS_PARAM));
        
        addOutput(createOutput<PJ301MPort>(Vec(153,329), module, FourierSynthesis::OUTPUT_LEFT));
        addOutput(createOutput<PJ301MPort>(Vec(182,329), module, FourierSynthesis::OUTPUT_RIGHT));

        // frequency display
        if (module) {
            auto* display = new FrequencyDisplay();
            display->box.pos = Vec(3 * RACK_GRID_WIDTH, 5 * RACK_GRID_WIDTH);
            display->box.size = Vec(9 * RACK_GRID_WIDTH, 3 * RACK_GRID_WIDTH);
            display->freqData = &module->freqMagnitudes;
            display->setNumBins(module->bufferSize / 2 + 1);
            addChild(display);
        }
    }
};

Model* modelFourierSynthesis = createModel<FourierSynthesis, FourierSynthesisWidget>("fourier-synthesis");
