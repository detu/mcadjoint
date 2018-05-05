//
// Created by Stefano Weidmann on 30.04.18.
//
#include "minimizer.hpp"
#include "sensitivity.hpp"
#include "vectorToBeMappedAsMatrix.hpp"
#include "logging.hpp"
#include <array>
#include "dumpToMatFile.hpp"
#include <stefCommonHeaders/xoroshiro.h>
#include <LBFGS.h>
#include <omp.h>

PermeabilitiesAndCost
matchWithPermeabilities(const FixedParameters& params, const int numberOfRows, const int numberOfCols,
                        const Real tolerance, const int maxIterations) {

    constexpr int seed = 42;

    const int numberOfCells = numberOfCols * numberOfRows;



    std::vector<Rng> rngs;

    for (int threadNumber = 0; threadNumber < omp_get_num_threads(); ++threadNumber) {
        rngs.emplace_back(seed, threadNumber);
    }



    const auto costFunctionForLbfgs = [&] (const Vector& logPermeabilitiesAsVector, Vector& sensitivities) -> Real {
        const Vector permeabilitiesAsVector = logPermeabilitiesAsVector.array().exp().matrix();
        const Eigen::Map<const Matrix> permeabilities(permeabilitiesAsVector.data(), numberOfRows, numberOfCols);

        const SensitivityAndCost sensitivityAndCost = computeSensitivityAndCost(params, permeabilities, rngs);

        dumpThis("sensitivities", sensitivityAndCost.sensitivity);
        dumpThis("permeabilities", permeabilities);
        dumpThis("cost", sensitivityAndCost.cost);

        writeToMatFile();

        sensitivities = sensitivityAndCost.sensitivity;

        return sensitivityAndCost.cost;
    };


    Vector logPermeabilitiesAsVector(numberOfCells);
    logPermeabilitiesAsVector.setZero();

    using namespace LBFGSpp;
    LBFGSParam<Real> lbfgsParam;
    lbfgsParam.epsilon = tolerance;
    lbfgsParam.max_iterations = maxIterations;
    lbfgsParam.min_step = 1e-9;

    LBFGSSolver<Real> lbfgsSolver(lbfgsParam);

    Real minimumCost = -1;


    (void) lbfgsSolver.minimize(costFunctionForLbfgs, logPermeabilitiesAsVector, minimumCost);
    PermeabilitiesAndCost permeabilitiesAndCost;

    const Vector permeabilitiesAsVector = logPermeabilitiesAsVector.array().exp().matrix();
    permeabilitiesAndCost.permeabilities = Eigen::Map<const Matrix>(permeabilitiesAsVector.data(), numberOfRows, numberOfCols);
    permeabilitiesAndCost.cost = minimumCost;

    dumpThis("permeabilities", permeabilitiesAndCost.permeabilities);
    dumpThis("cost", permeabilitiesAndCost.permeabilities);



    return permeabilitiesAndCost;




}
