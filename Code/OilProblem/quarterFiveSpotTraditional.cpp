//
// Created by Stefano Weidmann on 09.05.18.
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
#include "forward.hpp"

int n = -1;

spdlog::level::level_enum level = spdlog::level::level_enum::off;

void parseCommandLine(const int argc, const char** argv) {
    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    std::string levelName = "info";
    std::string matFileName = "fieldsQuarterFiveSpotTraditional.mat";
    cmdl({"-n", "--dimension"}) >> n;
    cmdl({"-l", "--level"}) >> levelName;
    cmdl({"-m", "--matfile"}) >> matFileName;

    level = spdlog::level::from_str(levelName);


    dumpInThisMatFile(matFileName);
}

int main(int argc, const char** argv) {
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

    params.finalTime = 0.1;
    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 1;
    };

    params.dynamicViscosityOil = 1;//0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 1;//0.0010518; // Water 20°C

    params.porosity = 1;
    params.initialSaturationsWater.resize(n, n);
    params.initialSaturationsWater.setConstant(0);

    const CellIndex drillCell = findDrillCell(n, n);
    drillCell(params.initialSaturationsWater) = 1;


    const Real milliDarcy = 1;
    params.initialPermeabilities.resizeLike(params.initialSaturationsWater);
    params.initialPermeabilities.setOnes();


    const int numberOfTimesteps = 2;

    const int stateSize = 2*n*n;
    Matrix adjointMatrix(stateSize * numberOfTimesteps, stateSize * numberOfTimesteps);
    Vector adjointRhs(stateSize * numberOfTimesteps);

    adjointRhs.setZero();
    adjointMatrix.setZero();

    SimulationState simulationState(n, n);

    for (int currentTimelevel = 0; currentTimelevel < numberOfTimesteps; ++currentTimelevel) {
        (void) stepForwardAndAdjointProblemTraditional(params, params.initialPermeabilities, currentTimelevel, simulationState, adjointMatrix, adjointRhs);
    }

    dumpThis("adjointMatrix", adjointMatrix);
    dumpThis("adjointRhs", adjointRhs);


    const Vector adjoint = adjointMatrix.lu().solve(adjointRhs);
    std::cout << "Adjoint =\n" << adjoint;
    dumpThis("adjoint", adjoint);
    writeToMatFile();

    return 0;

}