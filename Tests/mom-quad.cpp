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

double f(double x, double y)
{
    return sin(M_PI * M_PI * x * y);
}

int main()
{
    int n = 4;
    int q = 2 * n;

    BMoment2DQuad mom_quad(q, n);

    // sets quadrilateral to compute 
    double v1[2] = {0, 0};
    double v2[2] = {1, 0};
    double v3[2] = {1, 1};
    double v4[2] = {0, 1};
    mom_quad.setQuadrilateral(v1, v2, v3, v4);
    // you should set the quadrilateral before getting the integration points

    // sets function values
    arma::mat quadPoints = mom_quad.getIntegrationPoints();
    arma::vec Fval(quadPoints.n_rows, arma::fill::ones);

    for (int i = 0; i < Fval.n_rows; i++)
        Fval(i) = f(quadPoints(i, 0), quadPoints(i, 1));

    mom_quad.setFunction(Fval);

    // computes moments
    mom_quad.compute_moments();

    ofstream file;
    file.open("results/mom_quad2.txt");

    // prints moments
    file << "Bernstein moments computed with the function sin(pi * x)" << endl
         << endl;
    for (int i = 0; i < n + 1; i++)
        for (int j = 0; j < n + 1; j++)
            file << "Bmoment[" << i << ", " << j << "] : " << mom_quad.get_bmoment(i, j, 0) << endl;

    file.close();

    return 0;
}