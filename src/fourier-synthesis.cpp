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
    double* left_buffer_in;
    double* right_buffer_in;
    double* left_buffer_out;
    double* right_buffer_out;
    int bufferIndex;
    int sampleRateIndex;
    std::vector<double> leftFreqMagnitudes;
    std::vector<double> rightFreqMagnitudes;

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
        INPUT_SAMPLE_RATE,
        INPUT_WAVEFORM,
        INPUT_HARMONICS,
        NUM_INPUTS
    };

    enum OutputIds {
        OUTPUT_LEFT,
        OUTPUT_RIGHT,
        NUM_OUTPUTS
    };
    enum LightId {
        SIGNAL_LIGHT,
        NUM_LIGHTS
    };

    FourierSynthesis() {
        config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
        configParam(BUFFER_PARAM,1.f,12000.f,1.f,"Buffer Size");
		getParamQuantity(BUFFER_PARAM)->snapEnabled = true;
        configParam(SAMPLE_RATE_PARAM,0.f,50.f,1.f,"Sample rate reduction");
		// getParamQuantity(SAMPLE_RATE_PARAM)->snapEnabled = true;
        configParam(WAVEFORM_PARAM,0.f,2.f,1.f,"Waveform Type");
        configParam(HARMONICS_PARAM,1.f,100.f,10.f,"Number of Harmonics");
		getParamQuantity(HARMONICS_PARAM)->snapEnabled = true;
        // initialise params & buffers
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
        if (left_buffer_in) fftw_free(left_buffer_in);
        if (right_buffer_in) fftw_free(right_buffer_in);
        if (left_buffer_out) fftw_free(left_buffer_out);
        if (right_buffer_out) fftw_free(right_buffer_out);
    }

    void initialiseResources() {
        // free resources
        if (real_in) fftw_free(real_in);
        if (freq_out) fftw_free(freq_out);
        if (real_out) fftw_free(real_out);
        if (left_buffer_in) fftw_free(left_buffer_in);
        if (right_buffer_in) fftw_free(right_buffer_in);
        if (left_buffer_out) fftw_free(left_buffer_out);
        if (right_buffer_out) fftw_free(right_buffer_out);
        if (forward_plan) fftw_destroy_plan(forward_plan);
        if (backward_plan) fftw_destroy_plan(backward_plan);

        // allocate buffers 
        real_in = fftw_alloc_real(bufferSize);
        freq_out = fftw_alloc_complex(bufferSize / 2 + 1);
        real_out = fftw_alloc_real(bufferSize);
        left_buffer_in = fftw_alloc_real(bufferSize);
        right_buffer_in = fftw_alloc_real(bufferSize);
        left_buffer_out = fftw_alloc_real(bufferSize);
        right_buffer_out = fftw_alloc_real(bufferSize);

        // set real_out array to zero to prevent the module
        // from outputting undefined array content
        memset(real_out, 0, bufferSize * sizeof(double));
        memset(left_buffer_out, 0, bufferSize * sizeof(double));
        memset(right_buffer_out, 0, bufferSize * sizeof(double));

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

        if (inputs[INPUT_SAMPLE_RATE].isConnected()) sampleRate += inputs[INPUT_SAMPLE_RATE].getVoltage();
        if (inputs[INPUT_WAVEFORM].isConnected()) waveformType += inputs[INPUT_WAVEFORM].getVoltage();
        if (inputs[INPUT_HARMONICS].isConnected()) numHarmonics += inputs[INPUT_HARMONICS].getVoltage();

        if (sampleRateIndex < sampleRate) {
            sampleRateIndex++;
        } else {
            sampleRateIndex = 0;
            if (bufferIndex < bufferSize) {
                // stream data from the input to the buffer
                left_buffer_in[bufferIndex] = static_cast<double>(inputs[INPUT_LEFT].getVoltage());
                right_buffer_in[bufferIndex] = static_cast<double>(inputs[INPUT_RIGHT].getVoltage());
                // simultaneously output the data from the previous buffer
                // divide by bufferSize since the output of FFTW is scaled by the size of the input buffer

                outputs[OUTPUT_LEFT].setVoltage(left_buffer_out[bufferIndex] / bufferSize);
                outputs[OUTPUT_RIGHT].setVoltage(right_buffer_out[bufferIndex] / bufferSize);

                bufferIndex++;
            } else {
                // make sure you don't miss a sample here
                bufferIndex = 0;

                memcpy(real_in, left_buffer_in, bufferSize*sizeof(double));
                fftw_execute(forward_plan);
                // modify frequency domain with custom waveforms based on chosen wavetable
                applyCustomWaveform(waveformType, bufferSize, freq_out);
                // compute magnitudes for the frequency domain (upper display)
                leftFreqMagnitudes.resize(bufferSize / 2 + 1);
                for (int i = 0; i < bufferSize / 2 + 1; i++) {
                    leftFreqMagnitudes[i] = sqrt(freq_out[i][0] * freq_out[i][0] + freq_out[i][1] * freq_out[i][1]);
                }

                scale(leftFreqMagnitudes);
                fftw_execute(backward_plan);
                memcpy(left_buffer_out, real_out, bufferSize*sizeof(double));


                memcpy(real_in, right_buffer_in, bufferSize*sizeof(double));
                fftw_execute(forward_plan);
                // modify frequency domain with custom waveforms based on chosen wavetable
                applyCustomWaveform(waveformType, bufferSize, freq_out);
                // compute magnitudes for lower display
                rightFreqMagnitudes.resize(bufferSize / 2 + 1);
                for (int i = 0; i < bufferSize / 2 + 1; i++) {
                    rightFreqMagnitudes[i] = sqrt(freq_out[i][0] * freq_out[i][0] + freq_out[i][1] * freq_out[i][1]);
                }
                scale(rightFreqMagnitudes);
                fftw_execute(backward_plan);
                memcpy(right_buffer_out, real_out, bufferSize*sizeof(double));
                

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

        // copy DC component
        temp_freq_out[0][0] = freq_out[0][0]; 
        temp_freq_out[0][1] = freq_out[0][1];
        for (int bin = 1; bin < bufferSize / 2 + 1; ++bin) {
            if (waveformType != 0){
                // pre-emptive 90° rotation on base freq
                double temp = freq_out[bin][0];
                freq_out[bin][0] = -freq_out[bin][1];
                freq_out[bin][1] = temp;
            }
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
                if (waveformType <= 0) { // sine (no harmonics)
                    if (harmonic != 1) continue;
                    r = freq_out[bin][0];
                    i = freq_out[bin][1];
                } else {
                    if (waveformType <= 1) { // sawtooth
                        harmonicMagnitude /= harmonic;
                    } else if (waveformType > 1) { // square
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
                    r = harmonicCoefficient * harmonicMagnitude * cos(phase * harmonic + std::pow(-1,harmonic) * M_PI_2);
                    i = harmonicCoefficient * harmonicMagnitude * sin(phase * harmonic + std::pow(-1,harmonic) * M_PI_2);
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

    void scale(std::vector<double>& data) {
        double min_value = *std::min_element(data.begin(), data.end());
        double max_value = *std::max_element(data.begin(), data.end());

        // handle division by 0 error
        if (max_value - min_value == 0) return;

        for (auto& value : data) {
            if (std::isnan(value) || std::isinf(value)) {
                value = 0.0;  // Reset invalid values
            }
        }

        // normalize each value in the array
        for (auto& value : data) {
            value = (value - min_value) / (max_value - min_value);
        }
    }
};

struct FrequencyDisplay : TransparentWidget {
    FourierSynthesis* module;
    ModuleWidget* moduleWidget;
    std::vector<double>* leftFreqData = nullptr; // pointer to frequency data
    std::vector<double>* rightFreqData = nullptr;
    int numBins = 0;                         // number of bins in display

    void setNumBins(int bins) {
        numBins = bins;
    }

    void drawLayer(const DrawArgs& args, int layer) override {
        if (layer != 1)
            return; // Only draw on the foreground layer
        if (!leftFreqData || leftFreqData->empty()) return;
        if (!rightFreqData || leftFreqData->empty()) return;
        // get widget size
        NVGcontext* vg = args.vg;
        float width = box.size.x;
        float height = box.size.y;

        // update number of bins
        numBins = (*leftFreqData).size();
        (*leftFreqData).resize(numBins, 0.0f);
        (*rightFreqData).resize(numBins, 0.0f);

        // draw background (temporary for positioning)
        // nvgBeginPath(vg);
        // nvgRect(vg, 0, 0, width, height);
        // nvgFillColor(vg, nvgRGB(20, 20, 20));
        // nvgFill(vg);

        // normalise data
        // scale(*leftFreqData);
        // scale(*rightFreqData);
        // for (auto& value : *freqData) {
        //     value = std::log10(value);
        // }
        // nvgGlobalCompositeOperation(vg, NVG_LIGHTER);


        for (int pass = 0; pass < 3; pass++) {
            float alpha = (pass == 0) ? 0.2f : (pass == 1) ? 0.5f : 1.0f;
            float spread = (pass == 0) ? 6.0f : (pass == 1) ? 3.0f : 1.0f;
            // upper graph
            nvgBeginPath(vg);
            for (int i = 0; i < numBins; i++) {
                float x = (float)i / numBins * width;
                float barHeight = clamp((*leftFreqData)[i] * (height / 2), 0.0f, height / 2); // Scale frequency magnitude
                // nvgRect(vg, x, (height / 2) - barHeight, width / numBins, barHeight);
                nvgRect(vg, x - spread / 2, (height / 2) - barHeight - spread / 2, (width / numBins) + spread, barHeight + spread);
            }
            nvgFillColor(vg, nvgRGBA(230, 233, 169, alpha * 255));
            nvgFill(vg);

            // lower graph
            nvgBeginPath(vg);
            for (int i = 0; i < numBins; i++) {
                float x = (float)i / numBins * width;
                float barHeight = clamp((*rightFreqData)[i] * (height / 2), 0.0f, height / 2); // Scale frequency magnitude
                // nvgRect(vg, x, height / 2, width / numBins, barHeight);
                nvgRect(vg, x - spread / 2, (height / 2) - spread / 2, (width / numBins) + spread, barHeight + spread);
            }
            // nvgFillColor(vg, nvgRGB(230, 233, 169)); // LED white
            nvgFillColor(vg, nvgRGBA(230, 233, 169, alpha * 255));
            nvgFill(vg);
        }
        nvgGlobalCompositeOperation(vg, NVG_SOURCE_OVER); // Reset blending mode
        // nvgFillColor(vg, nvgRGB(0, 200, 255)); //cyan
        // nvgFillColor(vg, nvgRGB(0, 0, 0)); //black
        // nvgFillColor(vg, nvgRGB(230, 233, 169)); //LED white
        // nvgFillColor(vg, nvgRGB(255, 255, 255)); //white
    }
};


struct FourierSynthesisWidget : ModuleWidget {
    FrequencyDisplay* display;

    FourierSynthesisWidget(FourierSynthesis* module) {
        setModule(module);
        setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/fourier_bg.svg")));

        // addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
        // addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
        // addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
        // addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

        addInput(createInput<PJ301MPort>(Vec(18,329), module, FourierSynthesis::INPUT_LEFT));
        addInput(createInput<PJ301MPort>(Vec(47,329), module, FourierSynthesis::INPUT_RIGHT));
        addInput(createInput<PJ301MPort>(Vec(20,250), module, FourierSynthesis::INPUT_SAMPLE_RATE));
        addInput(createInput<PJ301MPort>(Vec(185,190), module, FourierSynthesis::INPUT_WAVEFORM));
        addInput(createInput<PJ301MPort>(Vec(120,250), module, FourierSynthesis::INPUT_HARMONICS));

        addParam(createParam<RoundLargeBlackKnob>(Vec(34,197), module, FourierSynthesis::BUFFER_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(57,235), module, FourierSynthesis::SAMPLE_RATE_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(135,197), module, FourierSynthesis::WAVEFORM_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(158,235), module, FourierSynthesis::HARMONICS_PARAM));
        
        addOutput(createOutput<PJ301MPort>(Vec(153,329), module, FourierSynthesis::OUTPUT_LEFT));
        addOutput(createOutput<PJ301MPort>(Vec(182,329), module, FourierSynthesis::OUTPUT_RIGHT));

        // addChild(createLight<MediumLight<RedLight>>(Vec(120, 290), module, FourierSynthesis::SIGNAL_LIGHT));

        // frequency display
        if (module) {
            auto* display = new FrequencyDisplay();
            display->box.pos = Vec(3 * RACK_GRID_WIDTH, 5 * RACK_GRID_WIDTH);
            display->box.size = Vec(9 * RACK_GRID_WIDTH, 3 * RACK_GRID_WIDTH);
            display->leftFreqData = &module->leftFreqMagnitudes;
            display->rightFreqData = &module->rightFreqMagnitudes;
            display->setNumBins(module->bufferSize / 2 + 1);
            display->module = module;
            display->moduleWidget = this;
            addChild(display);
        }
    }
};

Model* modelFourierSynthesis = createModel<FourierSynthesis, FourierSynthesisWidget>("fourier-synthesis");
