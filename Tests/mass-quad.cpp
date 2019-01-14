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

    // geometry element (default: [0,1]^2)
    Element<Element_t::QuadrilateralEl> el({{0, 0}, {1, 0}, {1, 1}, {0, 1}});

    BMass2DQuad mass_quad(q, n);

    // sets function values (1)
    {
        arma::vec Fval(q * q, arma::fill::ones);
        mass_quad.setFunctionValues(Fval);
    }

    // computes mass matrix
    mass_quad.computeMatrix();

    // prints mass matrix
    for (int i = 0; i < (n + 1) * (n + 1); i++)
    {
        for (int j = 0; j < (n + 1) * (n + 1); j++)
        {
            cout << scientific << mass_quad.getMatrixValue(i, j) << ", ";
        }
        cout << endl;
    }

    return 0;
}