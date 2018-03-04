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

using namespace std;

using namespace Eigen;

void Driver::solve_Burger() {
    const int n = 100;
    const int npar = n;
    double dx = 1.0 / double(n);
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

    Matrix<double, n, 1> b_loc, u0, u, uo, utarget;
    Matrix<double, n, n> A_dia, A_off, I_loc, P_off, w_off;
    Matrix<double, n, npar> c_loc;
    Matrix<double, npar, 1> E_D;

    b_loc *= 0.0;
    c_loc *= 0.0;
    I_loc *= 0.0;
    A_dia *= 0.0;
    A_off *= 0.0;
    P_off *= 0.0;
    w_off *= 0.0;
    E_D *= 0.0;

    int q = q_per_dof * n * npar;                      // number of random walks in total
    double* W;
    int* alpha_k;
    W = new double[q];
    alpha_k = new int[q];

    /******
     * targer final solution and initialization of initial solution solution
     ******/
    for (int i = 0; i < n; i++) {
        I_loc(i, i) = 1.0;
        u0(i) = double(i + 1) * double(n - i) * dx * dx;
        utarget(i) = u0(i) * 2.0;
    }
    /******
     * optimization loop
     ******/
    FILE* file;
    Gnuplot g = Gnuplot("gp");                                // setting up gnuplot process
    setTerminal(g);
    for (int k = 0; k < nopt; k++) {
        uo = u0;
        E_D *= 0.0;                                             // initialize estimator
        /******
         * time step loop
         ******/
        for (int jm = 0; jm < m; jm++) {                            // time step loop
            /******
             * HERE GOES DHE COMPUTATION OF THE DIAGONAL JACOBIAN BLOCK
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
                    if (i > 0) vl += u(i - 1);
                    double vr = u(i);
                    if (i < n - 1) vr += u(i + 1);
                    double aa, bb, cc;
                    if (vl > 0) aa = -1.0;
                    else aa = 0.0;
                    if (vr < 0) cc = 1.0;
                    else cc = 0.0;
                    if ((vl > 0) && (vr > 0)) bb = 1.0;
                    else if ((vl <= 0) && (vr <= 0)) bb = -1.0;
                    else bb = 0.0;
                    u(i) = uo(i) - dt / dx * bb * 0.5 * uo(i) * uo(i);
                    if (i > 0) u(i) -= dt / dx * aa * 0.5 * uo(i - 1) * uo(i - 1);
                    if (i < n - 1) u(i) -= dt / dx * cc * 0.5 * uo(i + 1) * uo(i + 1);
                    /******
                     * HERE GOES DHE COMPUTATION OF THE OFF-DIAGONAL JACOBIAN BLOCK
                     ******/
                    if (i > 0) A_off(i - 1, i) = aa * dt / dx * u(i - 1);
                    if (i < n - 1) A_off(i + 1, i) = cc * dt / dx * u(i + 1);
                    A_off(i, i) = bb * dt / dx * u(i) - 1.0;
                    /******
                     * HERE GOES DHE COMPUTATION OF THE B-VECTOR
                     ******/
                    if (jm == m - 2) {
                        b_loc(i) = 2.0 * (utarget(i) - u(i));
                    }
                    /******
                     * HERE GOES DHE COMPUTATION OF THE C-BLOCK
                     ******/
                    if (jm == 0) {
                        if (i > 0) c_loc(i, i - 1) = aa * dt / dx * u0(i - 1);
                        if (i < n - 1) c_loc(i, i + 1) = cc * dt / dx * u0(i + 1);
                        c_loc(i, i) = bb * dt / dx * u0(i) - 1.0;
                    }
                }
            }
            uo = u;
            /******
             * HERE GOES DHE COMPUTATION OF THE DIAGONAL JACOBIAN BLOCK
             ******/
            A_dia = I_loc;
            A_dia = I_loc - A_dia;
            A_off = -A_off;
            /******
             * do this except for the last time step
             ******/
            if (jm < m - 1) {
                /******
                 * transition probabilities and wights
                 ******/
                for (int i = 0; i < n; i++) {
                    double rowsum = 0.0;
                    for (int j = 0; j < n; j++) {
                        rowsum += abs(A_off(i, j));
                    }
                    for (int j = 0; j < n; j++) {
                        P_off(i, j) = abs(A_off(i, j)) / rowsum;              // local transition probabilities
                        if (A_off(i, j) > 0) w_off(i, j) = rowsum;           // local weights
                        else w_off(i, j) = -rowsum;           // ..
                    }
                }
            }
            /******
             * loop over random walks
             ******/
            for (int p = 0; p < q; p++) {                               // do the following for q random walks
                int alpha_k0 = p / (q_per_dof * npar);                      // start row index
                int jm0 = alpha_k0 / n;                              // first time step of random wal p
                int jpar = p % npar;                                  // parameter index
                if (jm >= jm0) {
                    if (jm ==
                        jm0) {                                        // do this only for the first step of random walk p
                        alpha_k[p] = alpha_k0;                                  // initial c componant of random walk p
                        W[p] = c_loc(alpha_k0 - jm0 * n, jpar) * double(n);            // initial W of random walk p
                        E_D[jpar] += W[p] * b_loc(alpha_k0 - jm0 * n);                // contribution to estimator
                    }
                    if (jm < m -
                             1) {                                         // do the following for time steps larger than jm0-1 and smaller than m-1
                        double r = rand.equal();                                // random colomn index
                        int alpha_kp1 = n * m;                                 // ..
                        double cum = 0.0;                                 // ..
                        for (int h = 0; h < n; h++) {                               // ..
                            cum += P_off(alpha_k[p] - jm * n, h);                      // ..
                            if (r < cum) {                                          // ..
                                alpha_kp1 = h + (jm + 1) * n;
                                break;                      // ..
                            }                                                     // ..
                        }
                        if (alpha_kp1 == n * m) break;                        // stop random walk, if absorbed
                        W[p] *= w_off(alpha_k[p] - jm * n, alpha_kp1 - (jm + 1) * n);  // update W
                        E_D[jpar] += W[p] * b_loc(alpha_kp1 - (jm + 1) * n);        // contribution to estimator
                        alpha_k[p] = alpha_kp1;                             // new row index becomes old colomn index
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
        E_D /= double(q / npar);                                    // average estimator
        u0 -= E_D * relax;
        /******
         * output for testing with gnuplot visualization
         ******/
        file = fopen("out", "w");
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
    delete W;
    cout << "done\n";
    pause();
    return;
}