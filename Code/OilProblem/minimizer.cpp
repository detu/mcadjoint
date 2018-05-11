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


PermeabilitiesAndCost
matchWithPermeabilities(const FixedParameters& params, const int numberOfRows, const int numberOfCols,
                        const Real tolerance, const int maxIterations) {

    constexpr int seed = 42;

    const int numberOfCells = numberOfCols * numberOfRows;

    Rng rng(seed);


    std::vector<Real> costHistory;
    std::vector<Matrix> permeabilitiesHistory;

    const auto costFunctionForLbfgs = [&] (const Vector& logPermeabilitiesAsVector, Vector& sensitivities) -> Real {
        log()->debug("logPermeabilitiesAsVector =\n{}", logPermeabilitiesAsVector);
        const Vector permeabilitiesAsVector = logPermeabilitiesAsVector.array().exp().matrix();
        const Eigen::Map<const Matrix> permeabilities(permeabilitiesAsVector.data(), numberOfRows, numberOfCols);

        const SensitivityAndCost sensitivityAndCost = computeSensitivityAndCost(params, permeabilities, rng);

        dumpThis("sensitivities", sensitivityAndCost.sensitivity);
        dumpThis("permeabilities", permeabilities);
        dumpThis("cost", sensitivityAndCost.cost);

        costHistory.push_back(sensitivityAndCost.cost);
        dumpThis("costHistory", costHistory);

        permeabilitiesHistory.push_back(permeabilities);
        dumpThis("permeabilitiesHistory", permeabilitiesHistory);


        writeToMatFile();

        sensitivities = sensitivityAndCost.sensitivity;

        return sensitivityAndCost.cost;
    };

    Matrix initialLogPermeabilities = params.initialPermeabilities.array().log().matrix();
    Vector logPermeabilitiesAsVector(Eigen::Map<Vector>(initialLogPermeabilities.data(), initialLogPermeabilities.size(), 1));


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
