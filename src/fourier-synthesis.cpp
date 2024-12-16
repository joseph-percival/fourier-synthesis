#include "fourier-synthesis-plugin.hpp"
#include <fftw3.h>
#include <vector>
#include <iostream>

struct FourierSynthesis : Module {
    int bufferSize;
    int sampleRate;
    int waveformType;
    int param3;
    fftw_plan forward_plan;
    // fftw_plan backward_plan;
    double* real_in;
    fftw_complex* freq_out;
    double* real_out;
    int bufferIndex;
    int sampleRateIndex;
    std::vector<double> freqMagnitudes;
    std::vector<std::vector<double>> wavetables;
    std::unordered_map<std::string, std::vector<std::vector<double>>> waveformCache;

    enum ParamIds {
        BUFFER_PARAM,
        SAMPLE_RATE_PARAM,
        WAVEFORM_PARAM,
        PARAM_3,
        NUM_PARAMS
    };

    enum InputIds {
        INPUT_LEFT,
        INPUT_RIGHT,
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
        configParam(BUFFER_PARAM,1.f,4000.f,1.f,"Buffer Size");
		getParamQuantity(BUFFER_PARAM)->snapEnabled = true;
        configParam(SAMPLE_RATE_PARAM,0.f,400.f,1.f,"Sample rate reduction");
		getParamQuantity(SAMPLE_RATE_PARAM)->snapEnabled = true;
        configParam(WAVEFORM_PARAM,0.f,2.f,1.f,"Waveform Type");
		getParamQuantity(WAVEFORM_PARAM)->snapEnabled = true;
		getParamQuantity(PARAM_3)->snapEnabled = true;
        //initialise params & buffers
        bufferSize = params[BUFFER_PARAM].getValue();
        sampleRate = params[SAMPLE_RATE_PARAM].getValue();
        waveformType = params[WAVEFORM_PARAM].getValue();
        param3 = params[PARAM_3].getValue();
        bufferIndex = 0;
        sampleRateIndex = 0;
        initialiseResources();
        // precompute wavetables
    }

        ~FourierSynthesis() {
        if (real_in) fftw_free(real_in);
        if (real_out) fftw_free(real_out);
        if (freq_out) fftw_free(freq_out);
        if (forward_plan) fftw_destroy_plan(forward_plan);
        // if (backward_plan) fftw_destroy_plan(backward_plan);
    }

    void initialiseResources() {
        // free resources
        if (real_in) fftw_free(real_in);
        if (freq_out) fftw_free(freq_out);
        if (real_out) fftw_free(real_out);
        if (forward_plan) fftw_destroy_plan(forward_plan);
        // if (backward_plan) fftw_destroy_plan(backward_plan);

        // allocate buffers
        real_in = fftw_alloc_real(bufferSize);
        freq_out = fftw_alloc_complex(bufferSize / 2 + 1);
        real_out = fftw_alloc_real(bufferSize);

        // set real_out array to zero to prevent the module
        // from outputting undefined array content
        memset(real_out, 0, bufferSize * sizeof(double));

        // generate plans
        forward_plan = fftw_plan_dft_r2c_1d(bufferSize, real_in, freq_out, FFTW_ESTIMATE);
        // backward_plan = fftw_plan_dft_c2r_1d(bufferSize, freq_out, real_out, FFTW_ESTIMATE);

        wavetables.clear();
        for (int i = 0; i < bufferSize / 2 + 1; ++i) {
            wavetables.push_back(generateWaveform(0, bufferSize, i));
        }
    }

    void process(const ProcessArgs& args) override {

        if (paramsModified()){
            //reevaluate params & buffers
            bufferSize = params[BUFFER_PARAM].getValue();
            sampleRate = params[SAMPLE_RATE_PARAM].getValue();
            waveformType = params[WAVEFORM_PARAM].getValue();
            param3 = params[PARAM_3].getValue();
            bufferIndex = 0;
            sampleRateIndex = 0;
            initialiseResources();
            precomputeWaveforms();
            // updateWavetables();
        }

        if (sampleRateIndex < sampleRate) {
            sampleRateIndex++;
        } else {
            sampleRateIndex = 0;
            if (bufferIndex < bufferSize) {
                // 
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
                // compute magnitudes for the frequency domain
                freqMagnitudes.resize(bufferSize / 2 + 1);
                for (int i = 0; i < bufferSize / 2 + 1; i++) {
                    freqMagnitudes[i] = sqrt(freq_out[i][0] * freq_out[i][0] + freq_out[i][1] * freq_out[i][1]);
                }


                // fftw_execute(backward_plan);
                // process data
                computeOutput();
            }
        }
    }

    bool paramsModified() {
        return params[BUFFER_PARAM].getValue() != bufferSize ||
               params[SAMPLE_RATE_PARAM].getValue() != sampleRate ||
               params[WAVEFORM_PARAM].getValue() != waveformType ||
               params[PARAM_3].getValue() != param3;
    }

    void computeOutput() {
        const auto& wavetables = waveformCache[getWaveformKey()];
        for (int i = 0; i < bufferSize; ++i) {
            real_out[i] = 0.0;
            for (int j = 0; j < bufferSize / 2 + 1; ++j) {
                double phase = std::atan2(freq_out[j][1], freq_out[j][0]);
                real_out[i] += freqMagnitudes[j] * wavetables[j][i] * std::cos(phase);
            }
        }
    }

    void precomputeWaveforms() {
        waveformCache.clear();
        for (int type = 0; type <= 2; ++type) {
            std::string key = getWaveformKey(type);
            waveformCache[key].resize(bufferSize / 2 + 1);
            for (int harmonic = 0; harmonic <= bufferSize / 2; ++harmonic) {
                waveformCache[key][harmonic] = generateWaveform(type, bufferSize, harmonic);
            }
        }
    }

    std::string getWaveformKey(int type = -1) const {
        if (type == -1) type = waveformType;
        return std::to_string(type) + "_" + std::to_string(bufferSize);
    }

    std::vector<double> generateWaveform(const int type, int size, int harmonic) {
        std::vector<double> table(size);
        for (int i = 0; i < size; ++i) {
            double phase = 2.0 * M_PI * i * harmonic / size;
            if (type == 0) // sine
                table[i] = sin(phase);
            else if (type == 1) // sawtooth
                table[i] = 2.0 * (phase / (2.0 * M_PI)) - 1.0; // normalized sawtooth
            else if (type == 2) // square
                table[i] = phase < M_PI ? 1.0 : -1.0;
            // maybe add more types?
        }
        return table;
    }
};

struct FrequencyDisplay : Widget {
    FourierSynthesis* module;
    std::vector<double>* freqData = nullptr; // pointer to frequency data
    int numBins = 0;                               // number of bins in display

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

        addParam(createParam<RoundLargeBlackKnob>(Vec(34,197), module, FourierSynthesis::BUFFER_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(57,235), module, FourierSynthesis::SAMPLE_RATE_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(135,197), module, FourierSynthesis::WAVEFORM_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(158,235), module, FourierSynthesis::PARAM_3));
        
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
