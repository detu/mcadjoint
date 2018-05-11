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



    spdlog::register_logger(stefCommonHeaders::setUpLog<stefCommonHeaders::NoMutex>(level));

    if (n < 0) {
        log()->error("Didn't specify a positive n. Are you sure you passed a positive value with the -n flag?");
        std::exit(1);
    }


    FixedParameters params;
    const Real atmosphericPressure = 1;
    params.overPressureDrill = [=] (const Real time) {
        return atmosphericPressure;
    };

    constexpr Real fieldWidth = 1;

    params.meshWidth = fieldWidth / Real(n);
    params.finalTime = 1;
    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 1;
    };

    params.dynamicViscosityOil = 1;//0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 1;//0.0010518; // Water 20°C

    params.porosity = 1;
    params.initialSaturationsWater.resize(n, n);
    params.initialSaturationsWater.setConstant(0);
    params.maxNumberOfTimesteps = 1e6;

    const CellIndex drillCell = findDrillCell(n, n);
    drillCell(params.initialSaturationsWater) = 1;


    const Real milliDarcy = 1;
    params.initialPermeabilities.resizeLike(params.initialSaturationsWater);

    constexpr bool useLognormal = false;
    constexpr bool constantOne = true;
    constexpr bool channel = false;
    if (useLognormal) {
        std::lognormal_distribution<Real> lognormalDistribution(milliDarcy, 1);
        Rng rng;

        for (int j = 0; j < params.initialPermeabilities.cols(); ++j) {
            for (int i = 0; i < params.initialPermeabilities.rows(); ++i) {
                params.initialPermeabilities(i, j) = lognormalDistribution(rng);
            }
        }
    } else if (constantOne) {
        params.initialPermeabilities.setOnes();
    } else if (channel) {
        params.initialPermeabilities.setConstant(1e-5);
        CellIndex pos = findDrillCell(n, n);
        const CellIndex last = findWellCell(n, n);

        do {
            pos(params.initialPermeabilities) = 1;
            pos.i--;
            pos.j++;
        } while (pos != last);

        last(params.initialPermeabilities) = 1;

    }




    (void) matchWithPermeabilities(params, n, n, 1e-6, 100);

    spdlog::drop_all();

    return 0;

}