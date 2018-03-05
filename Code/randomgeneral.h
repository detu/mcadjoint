/*
 *                        randomgeneral.h  -  description
 *                           -------------------
 *  begin                : July 2003
 *  copyright            : (C) 2003 by Patrick Jenny
 *  email                : jenny@ifd.mavt.ethz
 *
 *  this program is property of Patrick Jenny
 */

#ifndef RANDOMGENERAL_H
#define RANDOMGENERAL_H

#include <math.h>
#include <string>

using namespace std;

class RandomGeneral {
protected:
    bool   _iset;
    double _gset;
public:
    RandomGeneral() {
        _iset = false;
        _rmInv = 1.0 / double(RAND_MAX);
    }
    ~RandomGeneral() {};
    double equal() {
        return double(rand()) * _rmInv;
    }
    void gaussCorr(double corr, double& r1, double& r2) {
        double a = sqrt(1.0 + corr) * gauss();
        double b = sqrt(1.0 - corr) * gauss();
        r1       = cos(0.7853981635) * a - sin(0.7853981635) * b;
        r2       = sin(0.7853981635) * a + cos(0.7853981635) * b;
    }
    void gaussCorrLimit(double corr, double& r1, double& r2) {
        double a = sqrt(1.0 + corr) * gaussLimit();
        double b = sqrt(1.0 - corr) * gaussLimit();
        r1       = cos(0.7853981635) * a - sin(0.7853981635) * b;
        r2       = sin(0.7853981635) * a + cos(0.7853981635) * b;
    }
    double gauss() {
        double fac, v1, v2, r;
        if (!_iset) {
            _iset = true;
            do {
                v1 = 2.0 * equal() - 1.0;
                v2 = 2.0 * equal() - 1.0;
                r  = pow(v1, 2) + pow(v2, 2);
            } while ((r >= 1.0) || (r == 0.0));
            fac   = sqrt(-2.0 * log(r) / r);
            _gset  = v1 * fac;
            return  v2 * fac;
        } else {
            _iset = false;
            return _gset;
        }
    }
    double gaussLimit() {
        double g;
        do {
            g = gauss();
        } while (fabs(g) > 3);
        return g;
    }
    double gamma(double m, double v) {
        double alpha, beta, v1, v2, v3, x, y, s, g, e, argexp;
        if (m <  0.0)       {
            printf("Random::gamma: m<0");
            return 0.0;
        } else if (v <  0.0)       {
            printf("Random::gamma: v<0");
            return 0.0;
        } else if (v == 0.0)       {
            return m;
        } else if (v <  m * m * 1.e-8) {
            x = m + sqrt(v) * gaussLimit();
        } else {
            alpha = m * m / v;
            beta  = v / m;
            if (alpha < 1.0)       {
                printf("Random::gamma: alpha<0");
                return 0.0;
            } else {
                s = sqrt(2.0 * alpha - 1.0);
                do {
                    v3 = equal();
                    do {
                        do {
                            do {
                                v1 = 2.0 * equal() - 1.0;
                                v2 = 2.0 * equal() - 1.0;
                            } while ((v1 * v1 + v2 * v2 > 1.0) || (v1 == 0.0));
                            y = v2 / v1;
                            g = s * y + alpha - 1.0;
                        } while (g <= 0.0);
                        argexp = (alpha - 1.0) * log(g / (alpha - 1.0)) - s * y;
                    } while (argexp < -60.0);
                    e = (1.0 + y * y) * pow(2.718281828, argexp);
                } while (v3 > e);
                return beta * g;
            }
        }
        return 0.0;
    }
private:
    double _rmInv;
};

#endif

