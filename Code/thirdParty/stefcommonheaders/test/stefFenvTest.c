#include <stefCommonHeaders/stefFenv.h>
#include <stdio.h>
#include <math.h>

int main() {
    StefFenv_CrashOnFPEs(FE_ALL_EXCEPT & ~FE_INEXACT);


    double x = 0;
    printf("x =? ");
    scanf("%lf", &x);

    printf("exp(x) = %lf\n 1 / x = %lf, exp(-x) = %lf\n", exp(x), 1.0 / x, exp(-x));


    return 0;
}
