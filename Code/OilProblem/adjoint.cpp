//
// Created by Stefano Weidmann on 17.04.18.
//

#include "oilProblem.hpp"
#include <stefCommonHeaders/assert.h>
#include <random>
#include <stefCommonHeaders/xoroshiro.h>

[[deprecated]]
SparseMatrix makeMatrixColStochastic(SparseMatrix matrix) {
    for (int col = 0; col < matrix.cols(); ++col) {
        matrix.col(col) /= matrix.col(col).sum();
    }

    return matrix;
}


[[deprecated]]
Real getTransitionProbability(const int fromState, const int toState, const SparseMatrix& transitionProbabilitiesColStochastic) {
    return transitionProbabilitiesColStochastic.coeff(toState, fromState);
}


bool transitionState(RandomWalkState& currentState, const BVectorSurrogate& b,
                                const SparseMatrix& pressureResidualsByPressures,
                                const SparseMatrix& pressureResidualsBySaturationsWater,
                                const SparseMatrix& saturationsWaterResidualsByPressure,
                                const SparseMatrix& saturationsWaterResidualsBySaturationsWater,
                                const int numberOfRows, const int numberOfCols, Rng& rng) {

    bool stillInTheSameTimestep;
    if (currentState.cell == CellIndex::invalidCell()) {
        stillInTheSameTimestep = false;
        return stillInTheSameTimestep;
    }

    struct Candidate {
        const CellIndex cellIndex;
        const bool isAPressure;
        const Real correspondingEntryOfAMatrix;
        const Real correspondingEntryOfBVector;
    };

    std::vector<Real> candidateUnnormalizedProbabilities(10);
    std::vector<Candidate> candidates;
    candidates.reserve(10);

    constexpr bool iAmAPressure = true;
    constexpr bool iAmASaturation = !iAmAPressure;
    constexpr bool neighborIsAPressure = iAmAPressure;
    constexpr bool neighborIsASaturation = iAmASaturation;

    const SparseMatrix& pressureResidualsDerived = (currentState.isAPressure? pressureResidualsByPressures: pressureResidualsBySaturationsWater);
    const SparseMatrix& saturationWaterResidualsDerived = currentState.isAPressure? saturationsWaterResidualsByPressure: saturationsWaterResidualsBySaturationsWater;

    Real sumOfUnnormalizedProbabilities = 0;
    if (currentState.isAPressure) {
        // for pressures there's the probability to stay at the same place
        const CellIndex meToMyself = pressureToTransmissibilityIndex(currentState.cell, currentState.cell, numberOfRows);

        const Real correspondingEntryOfAMatrix = 1 - meToMyself(pressureResidualsByPressures);
        const Real correspondingEntryOfBVector = b(currentState.cell, iAmAPressure);
        const Real unnormalizedProbabilityOfStayingHere = std::abs(correspondingEntryOfAMatrix);
        ASSERT(unnormalizedProbabilityOfStayingHere > 0);

        Candidate candidate = {currentState.cell, iAmAPressure, correspondingEntryOfAMatrix, correspondingEntryOfBVector};
        candidates.push_back(std::move(candidate));
        candidateUnnormalizedProbabilities.push_back(unnormalizedProbabilityOfStayingHere);

        sumOfUnnormalizedProbabilities += unnormalizedProbabilityOfStayingHere;
    }

    // add neighbors
    for (const auto neighbor: currentState.cell.neighbors(numberOfRows, numberOfCols)) {
        const CellIndex neighborToMe = pressureToTransmissibilityIndex(neighbor, currentState.cell, numberOfRows);

        const Real correspondingEntryOfAMatrixPressure = -neighborToMe(pressureResidualsDerived);
        const Real correspondingEntryOfBVectorPressure = b(neighbor, neighborIsAPressure);
        const Real candidateUnnormalizedProbabilityPressure = std::abs(correspondingEntryOfAMatrixPressure);
        ASSERT(candidateUnnormalizedProbabilityPressure > 0);

        const Real correspondingEntryOfAMatrixSaturationWater = -neighborToMe(pressureResidualsDerived);
        const Real correspondingEntryOfBVectorSaturationWater = b(neighbor, neighborIsASaturation);
        const Real candidateUnnormalizedProbabilitySaturationWater = std::abs(correspondingEntryOfAMatrixSaturationWater);
        ASSERT(candidateUnnormalizedProbabilitySaturationWater > 0);

        Candidate candidatePressure = {neighbor, iAmAPressure, correspondingEntryOfAMatrixPressure, correspondingEntryOfBVectorPressure};
        candidates.push_back(std::move(candidatePressure));
        candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilityPressure);


        Candidate candidateSaturationWater = {neighbor, iAmASaturation, correspondingEntryOfAMatrixSaturationWater, correspondingEntryOfBVectorSaturationWater};
        candidates.push_back(std::move(candidateSaturationWater));
        candidateUnnormalizedProbabilities.push_back(candidateUnnormalizedProbabilitySaturationWater);


        sumOfUnnormalizedProbabilities += candidateUnnormalizedProbabilityPressure + candidateUnnormalizedProbabilitySaturationWater;

    }

    const Candidate absorptionCandidate = {CellIndex::invalidCell(), false, NAN, NAN };
    const Real unnormalizedAbsorptionProbability = std::abs(b(currentState.cell, currentState.isAPressure));

    if (unnormalizedAbsorptionProbability > 0) {
        candidates.push_back(absorptionCandidate);
        candidateUnnormalizedProbabilities.push_back(unnormalizedAbsorptionProbability);
    }


    std::discrete_distribution<int> choose(candidateUnnormalizedProbabilities.cbegin(), candidateUnnormalizedProbabilities.cend());
    const int nextIndex = choose(rng);
    const Candidate& chosenCandidate = candidates[nextIndex];

    if (chosenCandidate.cellIndex == absorptionCandidate.cellIndex) {
        stillInTheSameTimestep = false;
        return stillInTheSameTimestep;
    }

    stillInTheSameTimestep = true;
    if (!(currentState.isAPressure && chosenCandidate.isAPressure)) {
        ++currentState.currentTimelevel;
        stillInTheSameTimestep = false;
    }
    // update W (pg. 6199, top)
    currentState.W *= sumOfUnnormalizedProbabilities * chosenCandidate.correspondingEntryOfAMatrix;
    // update D (pg. 6199, top)
    currentState.D += currentState.W * chosenCandidate.correspondingEntryOfBVector;


    currentState.isAPressure = chosenCandidate.isAPressure;
    currentState.cell = chosenCandidate.cellIndex;
    return stillInTheSameTimestep;
}


std::vector<RandomWalkState> initializeRandomWalks(const FixedParameters& params, ConstMatrixRef permeabilities, SimulationState& initialSimulationState) {
    const Matrix totalMobilities = computeTotalMobilities(params.dynamicViscosityOil, params.dynamicViscosityWater, permeabilities, initialSimulationState.saturationsWater);
    const int numberOfCols = permeabilities.cols();
    const int numberOfRows = permeabilities.rows();
    const int numberOfParameters = permeabilities.size();
    const int numberOfCells = numberOfParameters;
    // solve pressure system
    Vector pressureRhs(numberOfCells);
    pressureRhs.setZero();

    const Real pressureDrillNow = params.overPressureDrill(initialSimulationState.time);
    const SparseMatrix pressureSystem = assemblePressureSystemWithBC(totalMobilities);
    const Real sourceAtDrillNow = std::abs(params.inflowPerUnitDepthWater(initialSimulationState.time));
    adaptRhsForPressure(sourceAtDrillNow, pressureRhs, initialSimulationState.saturationsWater.rows(), initialSimulationState.saturationsWater.cols());
    initialSimulationState.pressures = solvePressurePoissonProblem(pressureSystem, pressureRhs);

    const CellIndex drillCell = findDrillCell(numberOfRows, numberOfCols);
    const Real computedPressureAtDrillCell = drillCell(initialSimulationState.pressures);
    const Real measuredPressureAtDrillCell = params.overPressureDrill(0);

    const SparseMatrix pressureResidualsByPressures = computePressureResidualsDerivedByPressure(pressureSystem);


    const BVectorSurrogate b(computedPressureAtDrillCell, measuredPressureAtDrillCell, numberOfRows, numberOfCols);
    const Real firstTimestep = getFirstTimestep();
    const SparseMatrix pressureResidualsByLogPermeabilities = computePressureResidualsByLogPermeability(initialSimulationState.pressures, totalMobilities);

    const Matrix fluxFunctionFactors = computeFluxFunctionFactors(initialSimulationState.saturationsWater, params.porosity, params.dynamicViscosityWater, params.dynamicViscosityOil);

    const Matrix pressureDerivativesX = computeXDerivative(initialSimulationState.pressures, params.meshWidth);
    const Matrix darcyVelocitiesX = computeTotalDarcyVelocitiesX(totalMobilities, pressureDerivativesX);
    const Matrix fluxesX = computeFluxesX(fluxFunctionFactors, darcyVelocitiesX);

    const Matrix pressureDerivativesY = computeYDerivative(initialSimulationState.pressures, params.meshWidth);
    const Matrix darcyVelocitiesY = computeTotalDarcyVelocitiesY(totalMobilities, pressureDerivativesY);
    const Matrix fluxesY = computeFluxesY(fluxFunctionFactors, darcyVelocitiesY);

    const SparseMatrix saturationResidualsByLogPermeabilities = computeSaturationsWaterResidualsByLogPermeability(fluxesX, fluxesY, totalMobilities, firstTimestep, params.meshWidth);
    const CMatrixSurrogate c(pressureResidualsByLogPermeabilities, saturationResidualsByLogPermeabilities, numberOfRows, numberOfCols);

    return initializeRandomWalks(numberOfRows, numberOfCols, numberOfParameters, b, c);
}

std::vector<RandomWalkState> initializeRandomWalks(const int numberOfRows, const int numberOfCols, const int numberOfParameters, const BVectorSurrogate& b, const CMatrixSurrogate& c) {

    const int numberOfRandomWalksPerPressureCell = 5;

    const int numberOfRandomWalksPerParameter = numberOfRandomWalksPerPressureCell;
    const int numberOfRandomWalks = numberOfRandomWalksPerParameter * numberOfParameters;

    // We start randomWalksPerPressureCell random walks for every pressure state, which means
    // randomWalksPerPressureCell per cell
    // We hope that the saturation states are reached through the coupling of the equations
    
    std::vector<RandomWalkState> randomWalks(numberOfRandomWalks);
    for (int parameterIndex = 0; parameterIndex <  numberOfParameters; ++parameterIndex) {
        const CellIndex cell = CellIndex::fromLinearIndex(parameterIndex, numberOfRows);

        std::vector<CellIndex> neighbors;

        constexpr static std::array<CellIndex::Direction, 4> directionsToCheck = {
              CellIndex::Direction::EAST, CellIndex::Direction::WEST, CellIndex::Direction::NORTH, CellIndex::Direction::SOUTH
        };

        for (const auto direction: directionsToCheck) {
            if (cell.hasNeighbor(direction, numberOfRows, numberOfCols)) {
                neighbors.push_back(cell.neighbor(direction));
            }
        }

        for (int localRandomWalkIndex = 0; localRandomWalkIndex < numberOfRandomWalksPerPressureCell; ++localRandomWalkIndex) {
            const int whichCellIndex = localRandomWalkIndex % int(neighbors.size());
            const CellIndex neighborOrMyself = neighbors[whichCellIndex];

            constexpr bool wantAPressure = true;
            RandomWalkState initialState;
            initialState.cell = neighborOrMyself;
            initialState.isAPressure = wantAPressure;
            initialState.currentTimelevel = 0;
            initialState.W = numberOfRandomWalksPerPressureCell * c(neighborOrMyself, cell, wantAPressure);
            initialState.D = initialState.W * b(neighborOrMyself, wantAPressure);

            randomWalks.push_back(std::move(initialState));
        }
    }

    return randomWalks;



}
