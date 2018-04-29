//
// Created by Stefano Weidmann on 10.04.18.
//

#include "saturation.hpp"
#include "simulationState.hpp"
#include "logging.hpp"
#include "forward.hpp"
#include <EigenSimplematio.hpp>
#include "specialCells.hpp"
#include <signal.h>
#include <cstdlib>
#include <argh.h>
#include <stefCommonHeaders/logging.hpp>
#include <memory>


int n = -1;
std::unique_ptr<SimulationState> simulationState = nullptr;


void writeToFileAndExit(int signalNumber) {
    SMIO::EigenMatFile matFile("fieldsQuarterFiveSpot.mat");
    matFile.writeVariable("satWater", clamp(simulationState->saturationsWater, 0, 1));
    matFile.writeVariable("pressure", Matrix(simulationState->pressures.map));
    matFile.close();
    std::exit(signalNumber);
}


void parseCommandLine(const int argc, const char** argv) {
    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    std::string levelName;
    cmdl("-n") >> n;
    cmdl({"-l", "--level"}) >> levelName;

    LOGGER = stefCommonHeaders::setUpLog(spdlog::level::from_str(levelName));
    simulationState = std::unique_ptr<SimulationState>(new SimulationState(n, n));

}

int main(const int argc, const char** argv) {

    parseCommandLine(argc, argv);


    //feenableexcept(FE_ALL_EXCEPT & ~FE_INEXACT);


    const Real meshWidth = 1.0 / n;

    const CellIndex wellCell = findWellCell(n, n);
    const CellIndex drillCell = findDrillCell(n, n);

    const Matrix permeabilities = Matrix::Ones(n, n);

    FixedParameters params;

    params.dynamicViscosityOil = 0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 0.0010518; // Water 20°C
    params.finalTime = 0.1;

    const Real atmosphericPressure = 1e5;

    params.overPressureDrill = [=] (const Real time) {
        return 2*atmosphericPressure;
    };



    params.porosity = 1;


    simulationState->saturationsWater.setZero();
    params.meshWidth = meshWidth;


    params.inflowPerUnitDepthWater = [&] (const Real time) {
        return 3;
    };


    Vector pressureRhs(n*n);
    pressureRhs.setZero();

    signal(SIGINT, writeToFileAndExit);
    signal(SIGTERM, writeToFileAndExit);


    while (simulationState->time < params.finalTime) {
        const bool breakThroughHappened = stepForwardProblem(params, permeabilities, *simulationState);
        LOGGER->debug("saturations water =\n{}", simulationState->saturationsWater);
        LOGGER->info("time = {}", simulationState->time);
        if (breakThroughHappened) {
            LOGGER->info("Water broke though to well.");
            break;
        }
    }

    LOGGER->debug("dim sat water = ({}, {})", simulationState->saturationsWater.rows(), simulationState->saturationsWater.cols());

    writeToFileAndExit(0);
    return 0;
}