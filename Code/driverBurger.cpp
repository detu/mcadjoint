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
#include <cmath>

using namespace std;

using namespace Eigen;


// PROBLEM
// -------------------------------------------
// You are given a discretization of the Burgers equation
// \partial_t u + \partial_x (u^2/2) = 0
// u(x, 0) = u0(x)
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
    const int n = 100;
    const int npar = n;
    // dx: meshwidth for 1D Burger
    const double dx = 1.0 / double(n);
    // dt: timestep for 1D Burger
    const double dt = dx;

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

    VectorXd b_loc(VectorXd::Zero(n)), u0(VectorXd::Zero(n)), u(VectorXd::Zero(n)), uo(VectorXd::Zero(n)), utarget(
        VectorXd::Zero(n));
    MatrixXd A_dia(MatrixXd::Zero(n, n)), A_off(MatrixXd::Zero(n, n)), I_loc(MatrixXd::Zero(n, n)), P_off(
        MatrixXd::Zero(n, n)), w_off(n, n);
    MatrixXd c_loc(MatrixXd::Zero(n, npar));
    VectorXd E_D(VectorXd::Zero(npar));


    const int q = q_per_dof * n * npar;                      // number of random walks in total
    Eigen::VectorXd W(Eigen::VectorXd::Constant(q, 1, NAN));
    Eigen::VectorXi alpha_k(Eigen::VectorXi::Constant(q, 1, -1));

    /******
     * target final solution and initialization of initial solution solution
     ******/
    for (int i = 0; i < n; i++) {
        I_loc(i, i) = 1.0;
        u0(i) = double(i + 1) * double(n - i) * dx * dx;
        utarget(i) = u0(i) * 1.1;
    }
    /******
     * optimization loop
     ******/
    Gnuplot g = Gnuplot("gp");                                // setting up gnuplot process
    setTerminal(g);
    for (int k = 0; k < nopt; k++) {
        uo = u0;
        E_D *= 0.0;                                             // initialize estimator
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
            b_loc *= 0.0;
            c_loc *= 0.0;
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
                    /******
                     * HERE GOES THE COMPUTATION OF THE OFF-DIAGONAL JACOBIAN BLOCK
                     ******/
                    // Equations (34) on pg. 6198
                    if (i > 0) {
                        A_off(i - 1, i) = aa * dt / dx * u(i - 1);
                    }
                    A_off(i, i) = bb * dt / dx * u(i) - 1.0;
                    if (i < n - 1) {
                        A_off(i + 1, i) = cc * dt / dx * u(i + 1);
                    }
                    /******
                     * HERE GOES THE COMPUTATION OF THE B-VECTOR
                     ******/

                    // jm: timestep index
                    // m: number of timesteps
                    // b_loc = \partial F / \partial u
                    // F = sum over all k of (utarget(j, (m-2) * dt) - u(j, (m-2) * dt))^2
                    // b_loc(i) = \partial F / partial u_i
                    if (jm == m - 2) {
                        b_loc(i) = 2.0 * (utarget(i) - u(i));
                    }
                    /******
                     * HERE GOES THE COMPUTATION OF THE C-BLOCK
                     ******/
                    if (jm == 0) {
                        // c_loc is a matrix whose columns are the c vectors in c dot \psi used in the Monte Carlo approach.
                        // To compute the sensitivity, we want this matrix to be \partial Residual R / \partial parameters eta
                        // our parameters here are just the values of the initial solution u0 at the grid points

                        // we can reuse the formulas in equation (34) on pg. 6198 with t = 1
                        if (i > 0) {
                            c_loc(i, i - 1) = aa * dt / dx * u0(i - 1);
                        }
                        if (i < n - 1) {
                            c_loc(i, i + 1) = cc * dt / dx * u0(i + 1);
                        }
                        c_loc(i, i) = bb * dt / dx * u0(i) - 1.0;
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
                int alpha_k0 = p / (q_per_dof * npar);                      // start row index
                int jm0 = alpha_k0 / n;                              // first time step of random walk p. jm0 == p/q
                int jpar = p % npar;                                  // parameter index
                if (jm >= jm0) {
                    // this random walk has already started
                    if (jm == jm0) {
                        // this random walk has started just at this timestep
                        // do this only for the first step of random walk p
                        alpha_k[p] = alpha_k0;                                  // initial c component of random walk p
                        W[p] = c_loc(alpha_k0 - jm0 * n, jpar) * double(
                                   n);            // initial W of random walk p. Here the birth probability is 1/n for all states
                        E_D[jpar] += W[p] * b_loc(alpha_k0 - jm0 * n);                // contribution to estimator
                    }
                    if (jm < m - 1) {                                         // do the following for time steps larger than jm0-1 and smaller than m-1
                        double r = rand.equal();                                // random column index
                        int alpha_kp1 = n * m;                                 // ..
                        double cum = 0.0;                                 // ..
                        for (int h = 0; h < n; h++) {                               // ..
                            cum += P_off(alpha_k[p] - jm * n, h);                      // ..
                            if (r < cum) {                                          // ..
                                alpha_kp1 = h + (jm + 1) * n;
                                break;                      // ..
                            }                                                     // ..
                        }
                        if (alpha_kp1 == n * m) {
                            break;    // stop random walk, if absorbed
                        }
                        W[p] *= w_off(alpha_k[p] - jm * n, alpha_kp1 - (jm + 1) * n);  // update W
                        E_D[jpar] += W[p] * b_loc(alpha_kp1 - (jm + 1) * n);        // contribution to estimator
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
        // optimize by steepest descent
        u0 -= E_D * relax;
        cout << "||E_D|| = " << E_D.norm() << "\n";
        /******
         * output for testing with gnuplot visualization
         ******/
        FILE* file = fopen("out", "w");
        for (int i = 0; i < n; i++) {
            fprintf(file, "%le %le %le %le\n", double(i + 1) * double(n - i) * dx * dx, u0(i), u(i), utarget(i));
        }
        fclose(file);

        g.cmd("p[0:100][0:1]'out'u 1 w l,'out'u 2 w l,'out'u 3 w l,'out'u 4 w l");
        cout << "optimization step " << k << " completed\n";
    }
    /******
     * end of optimization loop
     ******/
    cout << "done\n";
    pause();
    return;
}