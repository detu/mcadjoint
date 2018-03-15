//
//  driver.cpp
//  Adjoint-MC
//
//  Created by Patrick Jenny on 14/02/18.
//  Copyright © 2018 Patrick Jenny. All rights reserved.
//

#include "driver.h"

#include <iostream>
#include <Eigen/Dense>
#include "randomgeneral.h"
#include "gnuplot_i.h"
#include "control.h"
#include "gnuPlotSetTerminal.hpp"
#include <stefCommonHeaders/assert.h>
#include <stefCommonHeaders/stefFenv.h>
#include <stefCommonHeaders/dbg.hpp>

#include <cmath>
#include <stdexcept>
#include <map>
#include <string>
#include <sstream>
#include "problemSelection.hpp"
#include "sparseDataSelection.hpp"

using namespace std;
using namespace Eigen;


// PROBLEM
// -------------------------------------------
// You are given a discretization of the Burgers equation
// \partial_t u + \partial_x (u^2/2) = 0
// u(x, 0) = u0(x)
// with zero Dirichlet boundary conditions at the boundaries 0, 1
// with forward euler and upwind.
// The corresponding residual is
// R_i^(t) := u_i^(t) - u_i^(t-1) + timestep / weshwidth * (f_(i+0.5)^(t-1) - f(i-0.5)^(t-1))
// where f_(i+0.5)^(t) = 0.5(u_i^(t))^2 if 0.5*(u_i^(t) + u_(i+1)^(t)) >= 0
//                       0.5(u_i+1^(t))^2 otherwise
// is the upwind flux.
// Given this discretization, find an initial vector u0_1, u0_n
// such that the solution closely matches a prescribed target solution utarget at a certain timestep T.
// Do this by minimizing
// F(u0) := 0.5 * sum over all cell indices i (u_i^(T) - utarget_i)^2


void Driver::solve_Burger() {
    StefFenv_CrashOnFPEs(FE_ALL_EXCEPT & ~FE_INEXACT);
    const int n = 100;
    // dx: meshwidth for 1D Burger
    const double dx = 1.0 / double(n);
    // dt: timestep for 1D Burger
    double dt = dx;

    RandomGeneral rand;
    Control control("control");
    double relax = 0.5;
    int m = 30;
    int q_per_dof = 5;                             // number of random walks per dof
    int nopt = 5;
    m = control.getInt("m");
    q_per_dof = control.getInt("q_per_dof");
    nopt = control.getInt("nopt");
    relax = control.getDouble("relax");

    const std::string problemDescription = control.getString("problem");
    const Problem problem = getProblemFromDescription(problemDescription);

    const std::string sparseDataDescription = control.getString("sparsedata");


    int npar;




    switch (problem) {
        case Problem::MATCH_FINAL_WITH_INITIAL:
        case Problem::MATCH_DATA_WITH_INITIAL: {
            npar = n;
            break;
        }

        case Problem::MATCH_DATA_WITH_VISCOSITY: {
            npar = 1;
            relax /= 10;
            dt /= 1000;
            break;
        }

        default: {
            throw logic_error("Not implemented!");
        }
    }


    cout << "Solving problem " << problemDescription << "\n";

    double viscosity = 0;
    VectorXd b_loc(VectorXd::Zero(n)), u(VectorXd::Zero(n)), uo(VectorXd::Zero(n));
    MatrixXd A_dia(MatrixXd::Zero(n, n)), A_off(MatrixXd::Zero(n, n)), I_loc(MatrixXd::Identity(n, n)), P_off(
        MatrixXd::Zero(n, n)), w_off(n, n);
    MatrixXd c_loc(MatrixXd::Zero(n, npar));
    VectorXd E_D(VectorXd::Zero(npar));


    const int q = q_per_dof * n * npar;                      // number of random walks in total
    VectorXd W(VectorXd::Constant(q, 1, NAN));
    VectorXi alpha_k(VectorXi::Constant(q, 1, -1));


    
    VectorXd u0(VectorXd::Zero(n));
    // First guess at initial solution
    for (int i = 0; i < n; i++) {
        u0(i) = double(i + 1) * double(n - i) * dx * dx;
    }

    Eigen::VectorXd utargetFinal;
    SparseData sparseData;
    switch (problem) {
        case Problem::MATCH_FINAL_WITH_INITIAL: {
            utargetFinal.resizeLike(u0);
            utargetFinal = u0 * 1.1;
            break;
        }

        case Problem::MATCH_DATA_WITH_VISCOSITY: {
            sparseData = getSparseDataFromDescription(sparseDataDescription);
            viscosity = 1;
            break;
        }
        case Problem::MATCH_DATA_WITH_INITIAL: {
            sparseData = getSparseDataFromDescription(sparseDataDescription);
            break;
        }

        default: {
            throw std::logic_error("Not implemented!");
        }
    }
    /******
     * optimization loop
     ******/
    Gnuplot g = Gnuplot("gp");                                // setting up gnuplot process
    setTerminal(g);

    double F;
    for (int k = 0; k < nopt; k++) {
        F = 0;
        uo = u0;
        E_D.setZero();                                           // initialize estimator
        /******
         * time step loop
         ******/
        // m: number of timesteps
        for (int jm = 0; jm < m; jm++) {                            // time step loop
            /******
             * HERE GOES THE COMPUTATION OF THE DIAGONAL JACOBIAN BLOCK
             *
             * A_dia = (\partial R(t)/\partial u(t))^T
             ******/
            b_loc.setZero();
            c_loc.setZero();
            for (int i = 0; i < n; i++) {
                /******
                 * compute solution u at new time level using upwind explicit Euler
                 ******/
                if (jm < m - 1) {

                    double vl = u(i);
                    if (i > 0) {
                        vl += u(i - 1);
                    }
                    double vr = u(i);
                    if (i < n - 1) {
                        vr += u(i + 1);
                    }

                    double aa, bb, cc; // the factors below equation (34) on pg. 6198

                    // a_i^(t)
                    if (vl > 0) {
                        aa = -1.0;
                    } else {
                        aa = 0.0;
                    }

                    // b_i^(t)
                    if ((vl > 0) && (vr > 0)) {
                        bb = 1.0;
                    } else if ((vl <= 0) && (vr <= 0)) {
                        bb = -1.0;
                    } else {
                        bb = 0.0;
                    }

                    // c_i^(t)
                    if (vr < 0) {
                        cc = 1.0;
                    } else {
                        cc = 0.0;
                    }

                    //  uo: u at old timelevel
                    u(i) = uo(i) - dt / dx * bb * 0.5 * uo(i) * uo(i);
                    if (i > 0) {
                        u(i) -= dt / dx * aa * 0.5 * uo(i - 1) * uo(i - 1);
                    }
                    if (i < n - 1) {
                        u(i) -= dt / dx * cc * 0.5 * uo(i + 1) * uo(i + 1);
                    }

                    double secondDerivative = 0;

                    if (problem == Problem::MATCH_DATA_WITH_VISCOSITY) {
                        // Boundary conditions are zero dirichlet

                        if (i == 0) {
                            secondDerivative = (u(1) - 2 * u(0)) / (dx * dx);
                        } else if (i == n-1) {
                            secondDerivative = (-2 * u(n-1) + u(n-2)) / (dx * dx);
                        } else {
                            secondDerivative = (u(i+1) - 2 * u(i) - u(i-1)) / (dx * dx);
                        }

                        const double viscousTerm = secondDerivative * viscosity;
                        u(i) += dt * viscousTerm;
                    }
                    /******
                     * HERE GOES THE COMPUTATION OF THE OFF-DIAGONAL JACOBIAN BLOCK
                     ******/
                    // Equations (34) on pg. 6198

                    if (i > 0) {
                        A_off(i - 1, i) = aa * dt / dx * u(i - 1);
                    }
                    //cout << u(i)<< "\n";
                    A_off(i, i) = bb * dt / dx * u(i) - 1.0;
                    if (i < n - 1) {
                        A_off(i + 1, i) = cc * dt / dx * u(i + 1);
                    }

                    if (problem == Problem::MATCH_DATA_WITH_VISCOSITY) {
                        const double viscousCorrection = viscosity * dt /(dx * dx);
                        if (i > 0) {
                            A_off(i-1, i) -= viscousCorrection;
                        }
                        A_off(i, i) += 2 * viscousCorrection;
                        if (i < n-1) {
                            A_off(i+1, i) -= viscousCorrection;
                        }
                    }


                    /******
                     * HERE GOES THE COMPUTATION OF THE B-VECTOR
                     ******/


                    // jm: timestep index
                    // m: number of timesteps
                    // b_loc = -\partial F / \partial u
                    // F = sum over all k of (utarget(j, (m-2) * dt) - u(j, (m-2) * dt))^2
                    // b_loc(i) = -\partial F / partial u_i

                    switch (problem) {
                        case Problem::MATCH_FINAL_WITH_INITIAL: {
                            if (jm == m - 2) {
                                b_loc(i) = 2.0 * (utargetFinal(i) - u(i));
                            } else {
                                b_loc(i) = 0;
                            }

                            break;
                        }

                        case Problem::MATCH_DATA_WITH_VISCOSITY:
                        case Problem::MATCH_DATA_WITH_INITIAL: {
                            switch (sparseData) {
                                case SparseData::ALWAYS_ZERO: {
                                    b_loc(i) = -2 * u(i);
                                    break;
                                }

                                case SparseData::HAT_PATTERN: {

                                    b_loc(i) = 0;
                                    if (2*i < n && (jm + 1 == i || jm + 1 == n - 1 - i)) {
                                        b_loc(i) = 2 * (1 - u(i));
                                    }
                                    break;
                                }

                                case SparseData::ZERO_DIAGONAL: {
                                    b_loc(i) = 0;
                                    if (i == jm) {
                                        b_loc(i) = -2 * u(i);
                                    }

                                    break;
                                }

                                case SparseData::LINEAR_OFFDIAGONAL: {
                                    b_loc(i) = 0;
                                    if (i == jm) {
                                        const double x = dx * i;
                                        b_loc(i) = 2 * (x - u(i));
                                    }

                                    break;
                                }

                                default: {
                                    throw std::logic_error("Not implemented!");
                                }
                            }

                            break;
                        }

                        default: {
                            throw std::logic_error("Not implemented!");
                        }
                    }

                    // Compute F value
                    switch (problem) {
                        case Problem::MATCH_FINAL_WITH_INITIAL: {
                            if (jm == m - 2) {
                                F += pow(utargetFinal(i) - u(i), 2);
                            }

                            break;
                        }

                        case Problem::MATCH_DATA_WITH_VISCOSITY:
                        case Problem::MATCH_DATA_WITH_INITIAL: {
                            switch (sparseData) {
                                case SparseData::ALWAYS_ZERO: {
                                    F += pow(u(i), 2);
                                    break;
                                }

                                case SparseData::HAT_PATTERN: {

                                    if (2*i < n && (jm + 1 == i || jm + 1 == n - 1 - i)) {
                                        F += std::pow(1 - u(i), 2);
                                    }
                                    break;
                                }

                                case SparseData::ZERO_DIAGONAL: {
                                    if (i == jm) {
                                        F += std::pow(u(i), 2);
                                    }
                                    break;
                                }

                                case SparseData::LINEAR_OFFDIAGONAL: {
                                    if (i < n - 1 && i == jm+1) {
                                        const double x = i * dx;
                                        F += std::pow(u(i) - x, 2);
                                    }
                                    break;
                                }

                                default: {
                                    throw std::logic_error("Not implemented!");
                                }
                            }

                            break;
                        }

                        default: {
                            throw std::logic_error("Not implemented!");
                        }
                    }

                    /******
                     * HERE GOES THE COMPUTATION OF THE C-BLOCK
                     ******/
                    if (jm == 0) {
                        // c_loc is a matrix whose columns are the c vectors in c dot \psi used in the Monte Carlo approach.
                        // To compute the sensitivity, we want this matrix to be \partial Residual R / \partial parameters eta
                        // our parameters here are just the values of the initial solution u0 at the grid points

                        // we can reuse the formulas in equation (34) on pg. 6198 with t = 1

                        switch (problem) {
                            case Problem::MATCH_DATA_WITH_INITIAL:
                            case Problem::MATCH_FINAL_WITH_INITIAL: {
                                if (i > 0) {
                                    c_loc(i, i - 1) = aa * dt / dx * u0(i - 1);
                                }
                                c_loc(i, i) = bb * dt / dx * u0(i) - 1.0;
                                if (i < n - 1) {
                                    c_loc(i, i + 1) = cc * dt / dx * u0(i + 1);
                                }
                                break;
                            }

                            case Problem::MATCH_DATA_WITH_VISCOSITY: {
                                DBG_PEXPR(secondDerivative);
                                DBG_PEXPR(dt);
                                DBG_PEXPR(i);
                                c_loc(i, 0) = -secondDerivative * dt;
                                break;
                            }

                            default: {
                                throw std::logic_error("Not implemented!");
                            }
                        }

                    }
                }
            }
            uo = u;
            /******
             * HERE GOES THE COMPUTATION OF THE DIAGONAL JACOBIAN BLOCK
             ******/

            A_dia.setZero();
            A_off = -A_off;
            /******
             * do this except for the last time step
             ******/
            if (jm < m - 1) {
                /******
                 * transition probabilities and weights
                 ******/
                // Here P = A with scaled rows such that they sum to 1.
                // ASK: Why not the optimal choice of equation (28) on pg. 6195?
                for (int i = 0; i < n; i++) {
                    double rowsum = 0.0;
                    for (int j = 0; j < n; j++) {
                        rowsum += abs(A_off(i, j));
                    }
                    for (int j = 0; j < n; j++) {
                        P_off(i, j) = abs(A_off(i, j)) / rowsum;              // local transition probabilities

                        // w_off is used to compute the capital W for our D estimator
                        // Compare to the formula on pg. 6190 and find that we implement it for our choice of P.
                        if (A_off(i, j) > 0) {
                            w_off(i, j) = rowsum;           // local weights.
                        } else {
                            w_off(i, j) = -rowsum;           // ..
                        }
                    }
                }
            }
            /******
             * loop over random walks
             ******/
            for (int p = 0; p < q; p++) {     // do the following for q random walks
                int alpha_k0 = p / (q_per_dof * npar);                    // start row index
                int jm0 = alpha_k0 / n;                              // first time step of random walk p. jm0 == p/q
                int jpar = p % npar;                                  // parameter index
                if (jm >= jm0) {
                    // this random walk has already started
                    if (jm == jm0) {
                        // this random walk has started just at this timestep
                        // do this only for the first step of random walk p
                        alpha_k[p] = alpha_k0;// initial c component of random walk p
                        DBG_PEXPR(jpar);
                        DBG_PEXPR(c_loc(alpha_k0, jpar));
                        W[p] = c_loc(alpha_k0, jpar) * double(n);            // initial W of random walk p. Here the birth probability is 1/n for all states
                        E_D[jpar] += W[p] * b_loc(alpha_k0);                // contribution to estimator
                        DBG_PEXPR(W[p]);
                        DBG_PEXPR(b_loc(alpha_k0));
                        DBG_PEXPR(alpha_k0);
                        DBG_PEXPR(p);
                        ASSERT(!isnan(E_D[jpar]));
                    }
                    if (jm < m - 1) {                                         // do the following for time steps larger than jm0-1 and smaller than m-1
                        double r = rand.equal();                                // random column index
                        int alpha_kp1 = n+1;                                 // ..
                        double cum = 0.0;                                 // ..
                        for (int h = 0; h < n; h++) {
                            // ..
                            cum += P_off(alpha_k[p], h);                      // ..
                            if (r < cum) {
                                alpha_kp1 = h;
                                break;                      // ..
                            }                                                     // ..
                        }
                        if (alpha_kp1 == n+1) {
                            break;    // stop random walk, if absorbed
                        }

                        W[p] *= w_off(alpha_k[p], alpha_kp1);  // update W
                        E_D[jpar] += W[p] * b_loc(alpha_kp1);        // contribution to estimator
                        ASSERT(!isnan(E_D[jpar]));

                        // advance state of
                        alpha_k[p] = alpha_kp1;                             // new row index becomes old column index
                    }
                }
            }
            /******
             * end of random walk loop
             ******/
        }
        /******
         * end of time step loop
         ******/
        E_D /= double(q / npar); // average estimator
        ASSERT(!isnan(E_D[0]));
        DBG_PEXPR(E_D(0));
        // optimize by steepest descent

         switch (problem) {
            case Problem::MATCH_FINAL_WITH_INITIAL:
            case Problem::MATCH_DATA_WITH_INITIAL: {
                u0 -= E_D * relax;
                break;
            }
            case Problem::MATCH_DATA_WITH_VISCOSITY: {
                viscosity -= E_D(0) * relax;
                viscosity = max(viscosity, 0.0);
                break;
            }

            default: {
                throw std::logic_error("Not implemented!");
            }
        }
        cout << "||E_D|| = " << E_D.norm() << "\nF = " << F <<  "\n";


        /******
         * output for testing with gnuplot visualization
         ******/
        FILE* file = fopen("out", "w");



        switch (problem) {
            case Problem::MATCH_FINAL_WITH_INITIAL: {
                for (int i = 0; i < n; i++) {
                    fprintf(file, "%le %le %le %le\n", double(i + 1) * double(n - i) * dx * dx, u0(i), u(i), utargetFinal(i));
                }

                fclose(file);

                g.cmd("p[0:100][0:1]'out'u 1 w l,'out'u 2 w l,'out'u 3 w l,'out'u 4 w l");
                break;
            }

            case Problem::MATCH_DATA_WITH_VISCOSITY: {
                for (int i = 0; i < n; i++) {
                    fprintf(file, "%le %le\n", double(i + 1) * double(n - i) * dx * dx, u(i));
                }

                fclose(file);

                g.cmd("p[0:100][0:1]'out'u 1 w l,'out'u 2 w l");
                std::cout << "viscosity = " << viscosity << "\n";
                break;
            }

            case Problem::MATCH_DATA_WITH_INITIAL: {

                for (int i = 0; i < n; i++) {
                    fprintf(file, "%le %le %le\n", double(i + 1) * double(n - i) * dx * dx, u0(i), u(i));
                }

                fclose(file);

                g.cmd("p[0:100][0:1]'out'u 1 w l,'out'u 2 w l, 'out'u 3 w l");
                break;
            }

            default: {
                throw std::logic_error("Not implemented!");
            }
        }

        cout << "optimization step " << k << " completed\n";
    }
    /******
     * end of optimization loop
     ******/
    cout << "done\n";
    pause();
    return;
}