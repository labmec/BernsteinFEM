/** mass-quad.cpp
 * 
 * This file makes a test for the computation
 * of the finite element Mass Matrix for
 * quadrilateral elements using Bernstein polynomials
 */

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <armadillo>

#include "MassM.h"
#include "JacobiGaussNodes.h"

using namespace std;

int main()
{
    int n = 4;
    int q = 2 * n;

    BMass2DQuad mass_quad(q, n);

    { // sets quadrilateral
        double v1[2] = {0, 0};
        double v2[2] = {1, 0};
        double v3[2] = {1, 1};
        double v4[2] = {0, 1};
        mass_quad.setQuadrilateral(v1, v2, v3, v4);
    }

    // sets function values (1)
    {
        arma::vec Fval(q * q, arma::fill::ones);
        mass_quad.setFunction(Fval);
    }

    // computes mass matrix
    mass_quad.compute_matrix();

    // prints moments computed for the mass matrix computation
    // for (int a1 = 0; a1 < n; a1++) {
    //     for (int a2 = 0; a2 < n - a1; a2 ++)
    //         cout << "Bmoment[" << a1 << ", " << a2 << ", " << n - a1 - a2 << "] : " << mass_quad.get_bmoment(a1, a2, 2 * n);
    // }

    // prints mass matrix
    for (int i = 0; i < (n+1) * (n+1); i++) {
        for (int j = 0; j < (n+1) * (n+1); j++) {
            cout << scientific << mass_quad.getMatrixValue(i, j) << ", ";
        }
        cout << endl;
    }

    return 0;
}