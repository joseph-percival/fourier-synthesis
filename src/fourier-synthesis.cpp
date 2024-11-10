#include "fourier-synthesis-plugin.hpp"
#include <iostream>
#include <fftw3.h>

struct FourierSynthesis : Module {
    int bufferSize;
    int numSinusoids;
    int param2;
    int param3;
    fftw_plan forward_plan;
    fftw_plan backward_plan;

    enum ParamIds {
        BUFFER_PARAM,
        SINES_PARAM,
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
        configParam(BUFFER_PARAM,0.f,400.f,1.f,"Buffer Size");
		getParamQuantity(BUFFER_PARAM)->snapEnabled = true;
        configParam(SINES_PARAM,0.f,40.f,1.f,"Number of Sinusoids");
    }

    void process(const ProcessArgs& args) override {
        // outputs[OUTPUT_LEFT].setVoltage(bufferSize);
        // outputs[OUTPUT_RIGHT].setVoltage(numSinusoids);
        if (paramsModified()){
            bufferSize = params[BUFFER_PARAM].getValue();
            numSinusoids = params[SINES_PARAM].getValue();
            param2 = params[PARAM_2].getValue();
            param3 = params[PARAM_3].getValue();

            int N = 8;
            std::vector<int> real_input = {10, 20, 30, 40, 50, 60, 70, 80};
            // Allocate input array for FFTW (type double for FFTW's precision)
            double* in = fftw_alloc_real(N);
            // Copy integer data to the FFTW input array as double
            for (int i = 0; i < N; i++) {
                in[i] = static_cast<double>(real_input[i]);
            }
            // Allocate output array for FFTW
            fftw_complex* out = fftw_alloc_complex(N / 2 + 1); // FFTW's real-to-complex output size
            generatePlan(N, in, out);
        }

        // Create a plan for real-to-complex 1D FFT
        // Execute the FFT
        fftw_execute(forward_plan);
        fftw_execute(backward_plan);
        // Print the output (frequency domain)
        // std::cout << "Frequency domain output:" << std::endl;
        // for (int i = 0; i < N; i++) {
        //     std::cout << in[i] << std::endl;
        // }
        std::cout << paramsModified();
        // // Need to scale output by N.
        // fftw_destroy_plan(forward_plan);
        // fftw_destroy_plan(backward_plan);
        // fftw_free(in);
        // fftw_free(out);
    }

    void generatePlan(int N, double* in, fftw_complex* out) {
        forward_plan = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
        backward_plan = fftw_plan_dft_c2r_1d(N, out, in, FFTW_ESTIMATE);
    }

    bool paramsModified() {
        if (params[BUFFER_PARAM].getValue() != bufferSize) return true;
        if (params[SINES_PARAM].getValue() != numSinusoids) return true;
        if (params[PARAM_2].getValue() != param2) return true;
        if (params[PARAM_3].getValue() != param3) return true;
        // else all paramaters are unmodified
        return false;
    }

};

struct FourierSynthesisWidget : ModuleWidget {
    FourierSynthesisWidget(FourierSynthesis* module) {
        setModule(module);
        setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/fourier_bg.svg")));

        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

        addInput(createInput<PJ301MPort>(Vec(18,329), module, FourierSynthesis::INPUT_LEFT));
        addInput(createInput<PJ301MPort>(Vec(47,329), module, FourierSynthesis::INPUT_RIGHT));

        addParam(createParam<RoundLargeBlackKnob>(Vec(34,197), module, FourierSynthesis::BUFFER_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(57,235), module, FourierSynthesis::SINES_PARAM));
        addParam(createParam<RoundLargeBlackKnob>(Vec(135,197), module, FourierSynthesis::PARAM_2));
        addParam(createParam<RoundLargeBlackKnob>(Vec(158,235), module, FourierSynthesis::PARAM_3));
        
        addOutput(createOutput<PJ301MPort>(Vec(153,329), module, FourierSynthesis::OUTPUT_LEFT));
        addOutput(createOutput<PJ301MPort>(Vec(182,329), module, FourierSynthesis::OUTPUT_RIGHT));
    }
};

Model* modelFourierSynthesis = createModel<FourierSynthesis, FourierSynthesisWidget>("fourier-synthesis");
