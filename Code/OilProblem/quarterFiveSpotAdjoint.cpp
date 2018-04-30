//
// Created by Stefano Weidmann on 30.04.18.
//
#include "minimizer.hpp"
#include "dumpToMatFile.hpp"
#include <argh.h>
#include "logging.hpp"
#include <stefCommonHeaders/logging.hpp>
#include <stefCommonHeaders/stefFenv.h>

int n = -1;

void parseCommandLine(const int argc, const char** argv) {
    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    std::string levelName = "info";
    std::string matFileName = "fieldsQuarterFiveSpotAdjoint.mat";
    cmdl("-n") >> n;
    cmdl({"-l", "--level"}) >> levelName;
    cmdl({"-m", "--matfile"}) >> matFileName;

    LOGGER = stefCommonHeaders::setUpLog(spdlog::level::from_str(levelName));

    if (n < 0) {
        LOGGER->error("Didn't specify a positive n. Are you sure you passed a positive value with the -n flag?");
        std::exit(1);
    }

    dumpInThisMatFile(matFileName);
}


int main(const int argc, const char** argv) {
    parseCommandLine(argc, argv);


    FixedParameters params;
    const Real atmosphericPressure = 1e5;
    params.overPressureDrill = [=] (const Real time) {
        return 2 * atmosphericPressure;
    };

    params.meshWidth = 200.0 / Real(n);
    params.finalTime = 1;
    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 3;
    };

    const Real densityOil = 920; // density of Heavy crude oil (https://wiki.anton-paar.com/en/crude-oil/)
    const Real kinematicViscosityOil = 2.86e-6; // Brent
    params.dynamicViscosityOil = densityOil * kinematicViscosityOil;
    params.dynamicViscosityWater = 0.0010518; // Water 20°C

    params.porosity = 1;
    params.initialSaturationsWater.resize(n, n);
    params.initialSaturationsWater.setZero();


    (void) matchWithPermeabilities(params, n, n, 1e-3, 100);

    return 0;

}