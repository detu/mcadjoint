//
// Created by Stefano Weidmann on 29.03.18.
//

#ifndef STEFCOMMONHEADERS_DARCYVELOCITYINTERNAL_HPP
#define STEFCOMMONHEADERS_DARCYVELOCITYINTERNAL_HPP

// Gradient zero at boundaries
__attribute__((pure))
static inline Matrix computeGradientComponent(ConstMatrixRef pressures, const Real meshWidth,
                                          const DerivativeDirection direction);

#endif //STEFCOMMONHEADERS_DARCYVELOCITYINTERNAL_HPP
