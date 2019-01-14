/** mom.cpp
 * 
 * This file makes a test for the computation
 * of Moments in 1D elements using the Bernstein
 * polynomials.
 * 
 * A Moment is the integral of a basis function
 * (Bernstein polynomials) multiplied by a weight function.
 */

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

    // geometry element (defaul [0, 1])
    // Element<Element_t::LinearEl> el({0, 1});

    // declares moments object 
    // BMoment1D mom1d(q, n, el);
    BMoment1D mom1d(q, n);

    auto points = mom1d.getIntegrationPoints();
    // sets function values
    {
        arma::vec Fval(points.n_rows, arma::fill::ones);

        for (int i = 0; i < q; i++)
            Fval(i) = f( (legendre_xi(q, i) + 1.0) * 0.5 );

        mom1d.setFunctionValues(Fval);
    }

    // computes moments
    mom1d.computeMoments();

    // prints moments
    cout << "Bernstein moments computed with the function sin(pi * x)" << endl
         << endl;
    for (int i = 0; i < n + 1; i++)
        cout << "Bmoment[" << i << "] : " << mom1d.getBMoment(i) << endl;

    return 0;
}