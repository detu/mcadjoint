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
#include <random>

int n = -1;

spdlog::level::level_enum level = spdlog::level::level_enum::off;

void parseCommandLine(const int argc, const char** argv) {
    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    std::string levelName = "info";
    std::string matFileName = "fieldsQuarterFiveSpotAdjoint.mat";
    cmdl({"-n", "--dimension"}) >> n;
    cmdl({"-l", "--level"}) >> levelName;
    cmdl({"-m", "--matfile"}) >> matFileName;

    level = spdlog::level::from_str(levelName);


    dumpInThisMatFile(matFileName);
}


int main(const int argc, const char** argv) {
    parseCommandLine(argc, argv);
    //StefFenv_CrashOnFPEs(FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW);

    auto sharedLogger = stefCommonHeaders::setUpLog(level);
    LOGGER = sharedLogger.get();

    if (n < 0) {
        LOGGER->error("Didn't specify a positive n. Are you sure you passed a positive value with the -n flag?");
        std::exit(1);
    }


    FixedParameters params;
    const Real atmosphericPressure = 1;
    params.overPressureDrill = [=] (const Real time) {
        return 2 * atmosphericPressure;
    };

    constexpr Real fieldWidth = 1;

    params.meshWidth = fieldWidth / Real(n);
    params.finalTime = 0.02;
    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 3;
    };

    params.dynamicViscosityOil = 0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 0.0010518; // Water 20°C

    params.porosity = 0.5;
    params.initialSaturationsWater.resize(n, n);
    params.initialSaturationsWater.setConstant(0);

    const CellIndex drillCell = findDrillCell(n, n);
    drillCell(params.initialSaturationsWater) = 1;


    const Real milliDarcy = 1;
    params.initialPermeabilities.resizeLike(params.initialSaturationsWater);

    std::lognormal_distribution<Real> lognormalDistribution(milliDarcy, 1);
    Rng rng;

    for (int j = 0; j < params.initialPermeabilities.cols(); ++j) {
        for (int i = 0; i < params.initialPermeabilities.rows(); ++i) {
            params.initialPermeabilities(i, j) = lognormalDistribution(rng);
        }
    }




    (void) matchWithPermeabilities(params, n, n, 1e-3, 100);

    return 0;

}