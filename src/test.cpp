#include "fourier-synthesis-plugin.hpp"

#include "../libs/fftw-3.3.10/config.h"

struct Test : Module {
    float channel[2] = {0.f,0.f};

    enum ParamIds {
        GAIN_PARAM,
        NUM_PARAMS
    };

    enum InputIds {
        ENUMS(INPUT,2),
        NUM_INPUTS
    };
    enum OutputIds {
        ENUMS(OUTPUT,2),
        NUM_OUTPUTS
    };
    enum LightId {
        NUM_LIGHTS
    };

    Test() {
        config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
        configParam(GAIN_PARAM,0.f,40.f,1.f,"Gain"); //first two are range, third is default
    }

    void process(const ProcessArgs& args) override {
        for (int i=0; i<2; i++){ // loop runs twice per sample
            if (inputs[INPUT+i].isConnected()){// does the input have a connection?
                channel[i] = inputs[INPUT+i].getVoltage() * 0.2 * params[GAIN_PARAM].getValue();
                channel[i] = channel[i] / pow(1 + pow(fabs(channel[i]), 2.5f), 0.4); // fabs = float absolute value
            } 
            else {
                channel[i] = 0.f;
            }
            outputs[OUTPUT + i].setVoltage(channel[i] * 5);
        }
    }
};
struct TestWidget : ModuleWidget {
    TestWidget(Test* module) {
        setModule(module);
        setPanel(APP->window->loadSvg(asset::plugin(pluginInstance, "res/vco1_panel.svg")));

        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 *RACK_GRID_WIDTH, 0)));
        addChild(createWidget<ScrewSilver>(Vec(RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));
        addChild(createWidget<ScrewSilver>(Vec(box.size.x - 2 * RACK_GRID_WIDTH, RACK_GRID_HEIGHT - RACK_GRID_WIDTH)));

        addInput(createInput<PJ301MPort>(Vec(11,350), module, Test::INPUT));
        addInput(createInput<PJ301MPort>(Vec(11+30,350), module, Test::INPUT+1));

        addParam(createParam<RoundHugeBlackKnob>(Vec(62,60), module, Test::GAIN_PARAM));
        
        addOutput(createOutput<PJ301MPort>(Vec(54,350), module, Test::OUTPUT));
        addOutput(createOutput<PJ301MPort>(Vec(54+30,350), module, Test::OUTPUT+1));
    }
};

Model* modelTest = createModel<Test, TestWidget>("test");
