#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <armadillo>
#include "Moments.h"
#include "JacobiGaussNodes.h"

using namespace std;

// testing if BMoments is functioning correctly

double f(double x)
{
    return sin(M_PI * x);
}

int main()
{
    int n = 4;
    int q = 2 * n;

    // interval variables
    double a = 0.0;
    double b = 1.0;

    BMoment2DTri mom_tri(q, n);

    // sets triangle to compute
    double v1[2] = {0, 0};
    double v2[2] = {1, 0};
    double v3[2] = {0, 1};
    mom_tri.setTriangle(v1, v2, v3);

    // sets function values
    {
        arma::vec Fval((q + 1) * (q + 1), arma::fill::ones); // it has ((q + 2) * (q + 1)) / 2 nodes, but is indexed using (q + 1) * (q + 1) positions

        // for (int i = 0; i < q * q; i++)
        //     Fval(i) = f( (legendre_xi(q, i) + 1.0) * 0.5 );

        mom_tri.setFunction(Fval);
    }

    // computes moments
    mom_tri.compute_moments();

    // prints moments
    cout << "Bernstein moments computed with the function sin(pi * x)" << endl
         << endl;
    for (int a1 = 0; a1 < (n + 1); a1++)
        for (int a2 = 0; a2 <= n - a1; a2++)
            cout << "Bmoment[" << a1 << ", " << a2 << ", " << n - a1 - a2 << "] : " 
            << mom_tri.get_bmoment(a1, a2, 0) /* this takes indexing int o account */
            << endl;

    // for (int i = 0; i < (n + 1) * (n + 1); i++)
    //     cout << "Bmoment[" << i << "] :\t" << mom_tri.get_bmoment(i) << endl;
    return 0;
}