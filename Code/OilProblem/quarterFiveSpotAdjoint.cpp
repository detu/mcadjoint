//
// Created by Stefano Weidmann on 30.04.18.
//
#include "minimizer.hpp"
#include "dumpToMatFile.hpp"
#include <argh.h>
#include "logging.hpp"
#include <stefCommonHeaders/logging.hpp>
#include <stefCommonHeaders/stefFenv.h>
#include "cellindex.hpp"
#include "specialCells.hpp"

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
    //StefFenv_CrashOnFPEs(FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW);


    FixedParameters params;
    const Real atmosphericPressure = 1;
    params.overPressureDrill = [=] (const Real time) {
        return 2 * atmosphericPressure;
    };

    constexpr Real fieldWidth = 1;

    params.meshWidth = fieldWidth / Real(n);
    params.finalTime = 1;
    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 3;
    };

    params.dynamicViscosityOil = 0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 0.0010518; // Water 20°C

    constexpr Real millidarcy = 1;
    params.porosity = 0.5;
    params.initialSaturationsWater.resize(n, n);
    params.initialSaturationsWater.setConstant(std::log(millidarcy));

    const CellIndex drillCell = findDrillCell(n, n);
    drillCell(params.initialSaturationsWater) = 1;



    (void) matchWithPermeabilities(params, n, n, 1e-3, 100);

    return 0;

}