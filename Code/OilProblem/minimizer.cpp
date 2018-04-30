//
// Created by Stefano Weidmann on 30.04.18.
//
#include "minimizer.hpp"
#include "sensitivity.hpp"
#include "vectorToBeMappedAsMatrix.hpp"
#include "logging.hpp"
#include <array>
#include "dumpToMatFile.hpp"

PermeabilitiesAndCost
matchWithPermeabilities(const FixedParameters& params, const int numberOfRows, const int numberOfCols,
                        const Real tolerance, const int maxIterations) {

    const int numberOfCells = numberOfCols * numberOfRows;

    std::array<VectorToBeMappedAsMatrix, 2> logPermeabilities = {
          VectorToBeMappedAsMatrix(Vector::Zero(numberOfCells), numberOfRows, numberOfCols),
          VectorToBeMappedAsMatrix(Vector::Zero(numberOfCells), numberOfRows, numberOfCols)
    };

    VectorToBeMappedAsMatrix& logPermeabilitiesCurrent = logPermeabilities[0];
    VectorToBeMappedAsMatrix& logPermeabilitiesOld = logPermeabilities[1];
    Matrix permeabilities(numberOfRows, numberOfCols);


    constexpr Real factorIfSuccessful = 1.1;
    constexpr Real factorIfUnsuccessful = 0.5;


    Real oldCost = INFINITY;
    Real lineSearchParameter = 1e-6;

    dumpThisOnExit("maxIterations", maxIterations);
    dumpThisOnExit("tolerance", tolerance);

    for (int iteration = 0; iteration < maxIterations; ++iteration) {
        LOGGER->info("iteration = {}", iteration);
        dumpThisOnExit("iteration", iteration);


        permeabilities = logPermeabilitiesCurrent.map.array().exp().matrix();
        const SensitivityAndCost sensitivityAndCost = computeSensitivityAndCost(params, permeabilities);

        LOGGER->debug("Sensitivities = {}", sensitivityAndCost.sensitivity);


        const Real normOfSensitivity = sensitivityAndCost.sensitivity.norm();
        dumpThisOnExit("converged", 0);
        dumpThisOnExit("sensitivities", sensitivityAndCost.sensitivity);

        if (normOfSensitivity < tolerance) {
            dumpThisOnExit("converged", 1);
            dumpThisOnExit("normOfSensitivity", normOfSensitivity);
            break;
        }

        const bool isFirstIteration = iteration == 0;
        if (sensitivityAndCost.cost >= oldCost) {
            // backtrack
            logPermeabilitiesCurrent = logPermeabilitiesOld.vec;
            lineSearchParameter *= factorIfSuccessful;
        } else {
            // advance
            oldCost = sensitivityAndCost.cost;
            logPermeabilitiesOld = logPermeabilitiesCurrent.vec;
            logPermeabilitiesCurrent.vec -= lineSearchParameter * sensitivityAndCost.sensitivity;
            lineSearchParameter *= factorIfSuccessful;
        }

        dumpThisOnExit("currentPermeabilities", permeabilities);
        dumpThisOnExit("cost", sensitivityAndCost.cost);


        LOGGER->info("cost = {}", sensitivityAndCost.cost);
        LOGGER->info("norm of sensitivity = {}", normOfSensitivity);
    }

    return {permeabilities, oldCost};


}
