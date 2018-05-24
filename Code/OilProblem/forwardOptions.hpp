//
// Created by Stefano Weidmann on 17.05.18.
//

#pragma once
#include "preconditioning.hpp"
const WhichPreconditioner preconditionerToUse = WhichPreconditioner::Q_TRANSPOSE_FROM_QR_AND_DIAGONAL;
constexpr bool startAddingRandomWalksAtBeginning = true;

