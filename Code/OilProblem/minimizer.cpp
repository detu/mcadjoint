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
#include <cppoptlib/problem.h>
#include <cppoptlib/solver/gradientdescentsolver.h>
#include <cppoptlib/solver/lbfgssolver.h>
#include <stefCommonHeaders/assert.h>


using namespace cppoptlib;

struct PermeabilitiesProblem: public Problem<Real> {
    using typename Problem<Real>::Scalar;
    using typename Problem<Real>::TVector;





    PermeabilitiesProblem(const FixedParameters& params, const int seed, const int maxIterations, const Real tolerance):
          params(params), maxIterations(maxIterations), tolerance(tolerance), rng(seed) {};


    Real value(const TVector& logPermeabilitiesAsRowVector) override {
        const Real foundCost = findSensitivityAndCost(logPermeabilitiesAsRowVector).cost;
        return foundCost;
    }


    SensitivityAndCost findSensitivityAndCost(const TVector& logPermeabilitiesAsRowVector) {
        const auto sensitivityAndCostIterator = sensitivitiesAndCosts.find(logPermeabilitiesAsRowVector);

        if (sensitivityAndCostIterator == sensitivitiesAndCosts.end()) {
            const int numberOfRows = params.initialPermeabilities.rows();
            const int numberOfCols = params.initialPermeabilities.cols();

            const Vector permeabilitiesAsVector = logPermeabilitiesAsRowVector.array().exp().matrix();
            const Eigen::Map<const Matrix> permeabilities(permeabilitiesAsVector.data(), numberOfRows, numberOfCols);
            const Eigen::Map<const Matrix> logPermeabilities(logPermeabilitiesAsRowVector.data(), numberOfRows, numberOfCols);

            const auto computedSensitivityAndCost = computeSensitivityAndCost(params, permeabilities, logPermeabilities,
                                                                              rng);

            sensitivitiesAndCosts.insert(std::make_pair(logPermeabilitiesAsRowVector, computedSensitivityAndCost));

            costHistory.push_back(computedSensitivityAndCost.cost);
            permeabilitiesHistory.push_back(permeabilities);


            dumpThis("costHistory", costHistory);
            dumpThis("permeabilitiesHistory", permeabilitiesHistory);
            dumpThis("sensitivity", computedSensitivityAndCost.sensitivity);
            dumpThis("cost", computedSensitivityAndCost.cost);
            dumpThis("permeabilities", permeabilities);
            writeToMatFile();



            return computedSensitivityAndCost;
        } else {
            return sensitivityAndCostIterator->second;
        }
    }

    void gradient(const TVector& logPermeabilitiesAsRowVector, TVector& sensitivity) override {
        sensitivity = findSensitivityAndCost(logPermeabilitiesAsRowVector).sensitivity.transpose();
    }


    // Hash function for Eigen matrix and vector.
    // The code is from `hash_combine` function of the Boost library. See
    // http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine .
    template<typename T>
    struct matrix_hash : std::unary_function<T, size_t> {
        std::size_t operator()(T const& matrix) const {
            // Note that it is oblivious to the storage order of Eigen matrix (column- or
            // row-major). It will give you the same hash value for two different matrices if they
            // are the transpose of each other in different storage order.
            size_t seed = 0;
            for (int i = 0; i < matrix.size(); ++i) {
                auto elem = *(matrix.data() + i);
                seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            }
            return seed;
        }
    };



    std::unordered_map<TVector, SensitivityAndCost, matrix_hash<TVector>> sensitivitiesAndCosts;
    const FixedParameters params;
    const int maxIterations;
    const Real tolerance;
    Rng rng;

    std::vector<Real> costHistory;
    std::vector<Matrix> permeabilitiesHistory;
};

PermeabilitiesAndCost
matchWithPermeabilities(const FixedParameters& params, const Real tolerance, const int maxIterations) {

    using namespace cppoptlib;
    PermeabilitiesProblem permeabilitiesProblem(params, 42, maxIterations, tolerance);



    const int numberOfCols = params.initialPermeabilities.cols();
    const int numberOfRows = params.initialPermeabilities.rows();

    GradientDescentSolver<PermeabilitiesProblem> solver;
    const Vector initialPermeabilitiesAsVector = Eigen::Map<const Vector>(params.initialPermeabilities.data(), params.initialPermeabilities.size());
    Vector logPermeabilitiesAsVector = initialPermeabilitiesAsVector.array().log().matrix();


    if (false) {
        log()->info(
              "Checking gradient implementation\n-----------------------------------------------------------------------");
        if (!permeabilitiesProblem.checkGradient(logPermeabilitiesAsVector, 1)) {
            throw std::logic_error("MC and FD gradients don't match!");
        }
        log()->info(
              "-----------------------------------------------------------------------\nGradient implementation OK");
    }

    solver.minimize(permeabilitiesProblem, logPermeabilitiesAsVector);

    const Real optimalCosts = permeabilitiesProblem.value(logPermeabilitiesAsVector);
    const Matrix optimalPermeabilities = Eigen::Map<Matrix>(logPermeabilitiesAsVector.data(), numberOfRows, numberOfCols).array().exp().matrix();

    dumpThis("optimalCosts", optimalCosts);
    dumpThis("optimalPermeabilities", optimalPermeabilities);
    writeToMatFile();

    return {optimalPermeabilities, optimalCosts};
}
