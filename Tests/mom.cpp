#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Moments.h"
#include "JacobiGaussNodes.h"

using namespace std;

// testing if BMoments is functioning correctly

double function(double x)
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

    BMoment1D mom1d(q, n);

    mom1d.setInterval(a, b);

    // sets function values
    {
        double *Fval = new double[q];
        for (int i = 0; i < q; i++)
            Fval[i] = 1.0; //function( (legendre_xi(q, i) + 1.0) * 0.5 );

        mom1d.setFunction(Fval);
        delete Fval;
    }

    // computes moments
    mom1d.compute_moments();

    // prints moments
    cout << "Bernstein moments computed with the function sin(pi * x)" << endl
         << endl;
    for (int i = 0; i < n + 1; i++)
        cout << "Bmoment[" << i << "] : " << mom1d.get_bmoment(i) << endl;

    return 0;
}