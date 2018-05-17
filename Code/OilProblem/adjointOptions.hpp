//
// Created by Stefano Weidmann on 17.05.18.
//

#pragma once
static constexpr bool alwaysAbsorbAtDrill = false;
static constexpr bool enableAbsorption = false;
static constexpr bool initializeJustAtBeginning = true;

static constexpr bool preferSaturations = false;
static constexpr Real preferenceForSaturations = 2;

static constexpr int numberOfRandomWalksToAdd = 1000;
static constexpr Real minimumProbabilityToBeAddedAtLeastOnce = 0.5 / Real(numberOfRandomWalksToAdd);