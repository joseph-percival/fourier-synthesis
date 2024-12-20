#include "fourier-synthesis-plugin.hpp"
#include <iostream>
#include <fftw3.h>

struct FourierSynthesis : Module {
    int bufferSize;
    int sampleRate;
    int param2;
    int param3;
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
        PARAM_2,
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
		getParamQuantity(PARAM_2)->snapEnabled = true;
		getParamQuantity(PARAM_3)->snapEnabled = true;
        //initialise params & buffers
        bufferSize = params[BUFFER_PARAM].getValue();
        sampleRate = params[SAMPLE_RATE_PARAM].getValue();
        param2 = params[PARAM_2].getValue();
        param3 = params[PARAM_3].getValue();
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

    void process(const ProcessArgs& args) override {
        // outputs[OUTPUT_LEFT].setVoltage(bufferSize);
        // outputs[OUTPUT_RIGHT].setVoltage(sampleRate);

        if (paramsModified()){
            crossfadeBuffer();
            //reevaluate params & buffers
            bufferSize = params[BUFFER_PARAM].getValue();
            sampleRate = params[SAMPLE_RATE_PARAM].getValue();
            param2 = params[PARAM_2].getValue();
            param3 = params[PARAM_3].getValue();
            bufferIndex = 0;
            sampleRateIndex = 0;
            initialiseResources();
        }

        if (sampleRateIndex < sampleRate) {
            sampleRateIndex++;
        } else {
            sampleRateIndex = 0;
            if (bufferIndex < bufferSize) {
                //stream data from the input to the buffer
                real_in[bufferIndex] = static_cast<double>(inputs[INPUT_LEFT].getVoltage());
                //simultaneously output the data from the previous buffer
                //divide by bufferSize since the output of FFTW is scaled by the size of the input buffer

                outputs[OUTPUT_LEFT].setVoltage(real_out[bufferIndex] / bufferSize);
                outputs[OUTPUT_RIGHT].setVoltage(inputs[INPUT_RIGHT].getVoltage());
                bufferIndex++;
            } else {
                //make sure you don't miss a sample here
                bufferIndex = 0;
                fftw_execute(forward_plan);
                // Compute magnitudes for the frequency domain
                freqMagnitudes.resize(bufferSize / 2 + 1);
                for (int i = 0; i < bufferSize / 2 + 1; i++) {
                    freqMagnitudes[i] = sqrt(freq_out[i][0] * freq_out[i][0] + freq_out[i][1] * freq_out[i][1]);
                }
                fftw_execute(backward_plan);
                //process data
            }
        }

        // std::cout << "Frequency domain output:" << std::endl;
        // for (int i = 0; i < N; i++) {
        //     std::cout << in[i] << std::endl;
        // }
        // std::cout << paramsModified();
    }

    void crossfadeBuffer() {
        for (int i = 0; i < bufferSize; i++) {
            real_out[i] = (real_out[i] * 0.5) + (real_in[i] * 0.5);
        }
    }

    bool paramsModified() {
        return params[BUFFER_PARAM].getValue() != bufferSize ||
               params[SAMPLE_RATE_PARAM].getValue() != sampleRate ||
               params[PARAM_2].getValue() != param2 ||
               params[PARAM_3].getValue() != param3;
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

};

#include <rack.hpp>
#include <vector>
using namespace rack;

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
        addParam(createParam<RoundLargeBlackKnob>(Vec(135,197), module, FourierSynthesis::PARAM_2));
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
