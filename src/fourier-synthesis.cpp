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
    double* real_in;
    fftw_complex* freq_out;
    double* real_out;
    int buffer_index;

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
        configParam(BUFFER_PARAM,1.f,400.f,1.f,"Buffer Size");
		getParamQuantity(BUFFER_PARAM)->snapEnabled = true;
        configParam(SINES_PARAM,0.f,40.f,1.f,"Number of Sinusoids");
        //initialise params & buffers
        bufferSize = params[BUFFER_PARAM].getValue();
        numSinusoids = params[SINES_PARAM].getValue();
        param2 = params[PARAM_2].getValue();
        param3 = params[PARAM_3].getValue();
        buffer_index = 0;
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
        // outputs[OUTPUT_RIGHT].setVoltage(numSinusoids);

        if (paramsModified()){
            //reevaluate params & buffers
            bufferSize = params[BUFFER_PARAM].getValue();
            numSinusoids = params[SINES_PARAM].getValue();
            param2 = params[PARAM_2].getValue();
            param3 = params[PARAM_3].getValue();
            buffer_index = 0;
            initialiseResources();
        }

        if (buffer_index < bufferSize) {
            //stream data from the input to the buffer
            real_in[buffer_index] = static_cast<double>(inputs[INPUT_LEFT].getVoltage());
            //simultaneously output the data from the previous buffer
            //divide by bufferSize since the output of FFTW is scaled by the size of the input buffer
            outputs[OUTPUT_LEFT].setVoltage(real_out[buffer_index] / bufferSize);
            outputs[OUTPUT_RIGHT].setVoltage(inputs[INPUT_RIGHT].getVoltage());
            buffer_index++;
        } else {
            //make sure you don't miss a sample here
            buffer_index = 0;
            fftw_execute(forward_plan);
            fftw_execute(backward_plan);
            //process data
        }

        // std::cout << "Frequency domain output:" << std::endl;
        // for (int i = 0; i < N; i++) {
        //     std::cout << in[i] << std::endl;
        // }
        // std::cout << paramsModified();
        // // Need to scale output by N.
    }

    bool paramsModified() {
        return params[BUFFER_PARAM].getValue() != bufferSize ||
               params[SINES_PARAM].getValue() != numSinusoids ||
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

        // generate plans
        forward_plan = fftw_plan_dft_r2c_1d(bufferSize, real_in, freq_out, FFTW_ESTIMATE);
        backward_plan = fftw_plan_dft_c2r_1d(bufferSize, freq_out, real_out, FFTW_ESTIMATE);
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
