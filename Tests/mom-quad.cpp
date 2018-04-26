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
// you can change the function to use in the computations by
// modifying the function below, and uncommenting the code under
// the "sets function values" section

double f(double x)
{
    return sin(M_PI * x);
}

int main()
{
    int n = 4;
    int q = 2 * n;

    BMoment2DQuad mom_quad(q, n);

    // sets triangle to compute
    double v1[2] = {0, 0};
    double v2[2] = {1, 0};
    double v3[2] = {1, 1};
    double v4[2] = {0, 1};
    mom_quad.setQuadrilateral(v1, v2, v3, v4);

    // sets function values
    arma::vec Fval(q * q, arma::fill::ones);

    // for (int i = 0; i < q * q; i++)
    //     Fval(i) = f( (legendre_xi(q, i) + 1.0) * 0.5 );

    mom_quad.setFunction(Fval);

    // computes moments
    mom_quad.compute_moments();

    // prints moments
    cout << "Bernstein moments computed with the function sin(pi * x)" << endl
         << endl;
    for (int i = 0; i < n + 1; i++)
        for (int j = 0; j < n + 1; j++)
            cout << "Bmoment[" << i << ", " << j << "] : " << mom_quad.get_bmoment(i, j, 0) << endl;

    return 0;
}