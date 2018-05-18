//
// Created by Stefano Weidmann on 09.05.18.
//

#ifndef JUST_COMPUTE_ADJOINT
    #error "MUST DEFINE JUST_COMPUTE_ADJOINT"
#endif

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
#include "utils.hpp"
#include "sensitivity.hpp"
#include <vector>

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

    const int numberOfTimesteps = 20;
    params.porosity = 1;
    params.initialSaturationsWater.resize(n, n);
    params.initialSaturationsWater.setConstant(0);
    params.maxNumberOfTimesteps = numberOfTimesteps;

    const CellIndex drillCell = findDrillCell(n, n);
    drillCell(params.initialSaturationsWater) = 1;


    const Real milliDarcy = 1;
    params.initialPermeabilities.resizeLike(params.initialSaturationsWater);
    params.initialPermeabilities.setOnes();



    const int stateSize = 2*n*n;


    adjointRhs.setZero();
    adjointMatrix.setZero();


    SimulationState simulationState(n, n);
    simulationState.saturationsWater = params.initialSaturationsWater;
    SimulationState mcSimulationState(n, n);


    Rng rng(88);
    const auto sensitivityAndCost = computeSensitivityAndCost(params, params.initialPermeabilities, params.initialPermeabilities.array().log().matrix(),
                                                              rng);

    dumpThis("adjointMC", sensitivityAndCost.sensitivity);

    ASSERT(allFinite(adjointMatrix));
    ASSERT(allFinite(adjointRhs));

    dumpThis("adjointMatrixTrad", adjointMatrix);
    dumpThis("adjointMatrixTradEigenvalues", Vector(adjointMatrix.eigenvalues().array().abs().matrix()));
    dumpThis("adjointRhsTrad", adjointRhs);
    //dumpThis("correctedAdjointMatrixTradEigenvalues", (adjointMatrix.diagonal().asDiagonal().inverse() * adjointMatrix).eigenvalues().array().abs().matrix());


    const Vector adjoint = adjointMatrix.lu().solve(adjointRhs);
    dumpThis("adjointTrad", adjoint);


    writeToMatFile();

    return 0;

}