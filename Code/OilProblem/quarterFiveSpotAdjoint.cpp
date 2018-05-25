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
#include "sensitivity.hpp"
#include <random>

int n = -1;

spdlog::level::level_enum level = spdlog::level::level_enum::off;
int maxNumberOfTimesteps = -1;
bool compareSensitivities;
int numberOfRandomWalksToAdd = -1;

bool enableAntitheticSampling = false;
Real regularizationParameter = 0;
std::string logFileName = "log.log";
bool symmetrizeGradient = false;
bool traditionalMinimization = false;

void parseCommandLine(const int argc, const char** argv) {
    argh::parser cmdl;
    cmdl.parse(argc, argv, argh::parser::PREFER_PARAM_FOR_UNREG_OPTION);

    std::string levelName = "info";
    std::string matFileName = "fieldsQuarterFiveSpotAdjoint.mat";

    compareSensitivities = bool(cmdl[{"-s"}]);
    enableAntitheticSampling = bool(cmdl[{"-a", "--antithetic-sampling"}]);
    symmetrizeGradient = bool(cmdl[{"-S", "-symmetrize-gradient"}]);
    traditionalMinimization = bool(cmdl[{"-T", "--traditional"}]);
    cmdl({"-r", "--random-walks"}) >> numberOfRandomWalksToAdd;
    cmdl({"-n", "--dimension"}) >> n;
    cmdl({"-l", "--level"}) >> levelName;
    cmdl({"-m", "--matfile"}) >> matFileName;
    cmdl({"-M", "--max-timesteps"}) >> maxNumberOfTimesteps;
    cmdl({"-p", "--regularization-parameter"}) >> regularizationParameter;
    cmdl({"-L", "--log"}) >> logFileName;

    level = spdlog::level::from_str(levelName);


    dumpInThisMatFile(matFileName);
}


int main(const int argc, const char** argv) {
    parseCommandLine(argc, argv);
    //StefFenv_CrashOnFPEs(FE_ALL_EXCEPT & ~FE_INEXACT & ~FE_UNDERFLOW);
    SMIO_DisableCompression();



    spdlog::register_logger(stefCommonHeaders::setUpLog<stefCommonHeaders::NoMutex>(level, true, logFileName.c_str()));

    if (n <= 0) {
        log()->error("Didn't specify a positive n. Are you sure you passed a positive value with the -n flag?");
        std::exit(1);
    }

    if (maxNumberOfTimesteps <= 0) {
        log()->error("Didn't specify a positive maximum number of timesteps. Are you sure you passed a positive value with the -M flag?");
        std::exit(1);
    }

    if (numberOfRandomWalksToAdd <= 0) {
        log()->error("Didn't specify a positive maximum number of random walks to add. Are you sure you passed a positive value with the -r flag?");
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

    params.regularizationParameter = regularizationParameter;

    params.dynamicViscosityOil = 0.630; // SAE motor oil 20°C
    params.dynamicViscosityWater = 0.0010518; // Water 20°C

    params.porosity = 1;
    params.initialSaturationsWater.resize(n, n);
    params.initialSaturationsWater.setConstant(0);
    findDrillCell(n, n)(params.initialSaturationsWater) = 1;
    params.maxNumberOfTimesteps = maxNumberOfTimesteps;
    params.numberOfRandomWalksToAdd = numberOfRandomWalksToAdd;
    params.enableAntitheticSampling = enableAntitheticSampling;
    params.symmetrizeGradient = symmetrizeGradient;
    params.traditionalMinimization = traditionalMinimization;


    if (params.traditionalMinimization) {
        log()->info("Traditional minimization");
    } else {
        log()->info("MC based minimization");
    }

    const Real milliDarcy = 1;
    params.initialPermeabilities.resizeLike(params.initialSaturationsWater);

    constexpr bool useLognormal = false;
    constexpr bool constant = true;
    constexpr bool channel = false;
    if (useLognormal) {
        std::lognormal_distribution<Real> lognormalDistribution(milliDarcy, 1);
        Rng rng;

        for (int j = 0; j < params.initialPermeabilities.cols(); ++j) {
            for (int i = 0; i < params.initialPermeabilities.rows(); ++i) {
                params.initialPermeabilities(i, j) = lognormalDistribution(rng);
            }
        }
    } else if (constant) {
        params.initialPermeabilities.setConstant(1);
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


    if (compareSensitivities) {
        log()->info("Just comparing sensitivities");

        const Matrix permeabilitiesToTest = params.initialPermeabilities;
        const Matrix logPermeabilitiesToTest = permeabilitiesToTest.array().log().matrix();
        Rng compareSensitivitiesRng;

        const SensitivityAndCost sensitivityAndCostMC = computeSensitivityAndCost(params,
                                                                                  permeabilitiesToTest, logPermeabilitiesToTest,
                                                                                  compareSensitivitiesRng);

        const SensitivityAndCost sensitivityAndCostTrad = computeSensitivityAndCostTraditional(params, permeabilitiesToTest, logPermeabilitiesToTest);

        const Vector sensitivityDiscrepancies = (sensitivityAndCostMC.sensitivity - sensitivityAndCostTrad.sensitivity).cwiseAbs();
        int sensitivityMaxAbsDiffLoc = -1;
        const Real sensitivityMaxAbsDiff = sensitivityDiscrepancies.maxCoeff(&sensitivityMaxAbsDiffLoc);
        const Real sensitivityAvgAbsDiff = sensitivityDiscrepancies.mean();
        log()->info("Cost traditional = {}, cost MC = {} (should be equal!)", sensitivityAndCostTrad.cost, sensitivityAndCostMC.cost);

        log()->info("Sensitivities traditional =\n{}\n\nSensitivities MC =\n{}\n\nSensitivity avg abs diff = {}\nSensitivity max abs diff = {} at index {} between {} for traditional and {} for MC",
                    sensitivityAndCostTrad.sensitivity, sensitivityAndCostMC.sensitivity, sensitivityAvgAbsDiff, sensitivityMaxAbsDiff, sensitivityMaxAbsDiffLoc,
                    sensitivityAndCostTrad.sensitivity(sensitivityMaxAbsDiffLoc), sensitivityAndCostMC.sensitivity(sensitivityMaxAbsDiffLoc));


        dumpThis("sensitivityDiscrepancies", sensitivityDiscrepancies);
        dumpThis("sensitivityMaxAbsDiffLoc", Real(sensitivityMaxAbsDiffLoc));
        dumpThis("sensitivityMaxAbsDiff", sensitivityMaxAbsDiff);
        dumpThis("sensitivityAvgAbsDiff", sensitivityAvgAbsDiff);
        writeToMatFile();
    } else {
        (void) matchWithPermeabilities(params, 1e-6, 100);

    }


    spdlog::drop_all();

    return 0;

}