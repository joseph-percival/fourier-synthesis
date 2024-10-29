#include "fourier-synthesis-plugin.hpp"

#include "../libs/fftw-3.3.10/config.h"

struct Test : Module {
    enum ParamIds {
        NUM_PARAMS
    };

    enum InputIds {
        ENUMS(INPUT,2),
        NUM_INPUTS
    };
    enum OutputIds {
        ENUMS(OUTPUT,1),
        NUM_OUTPUTS
    };
    enum LightId {
        NUM_LIGHTS
    };
    Test() {
        config(NUM_PARAMS, NUM_INPUTS, NUM_OUTPUTS, NUM_LIGHTS);
    }

    void process(const ProcessArgs& args) override {

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
        
        addOutput(createOutput<PJ301MPort>(Vec(54,350), module, Test::OUTPUT));
    }
};

Model* modelTest = createModel<Test, TestWidget>("test");
