/** mas.cpp
 * 
 * This file makes a test for the computation
 * of the finite element Mass Matrix for 1D elements
 * using Bernstein polynomials
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

    // geometry element (default: [0,1])
    // Element<Element_t::LinearEl> el({0, 1});

    // Mass matrix object
    // BMass1D mass1d(q, n, el);
    BMass1D mass1d(q, n);

    // sets function values (1)
    {
        arma::vec Fval(q, arma::fill::ones);
        mass1d.setFunctionValues(Fval);
    }

    // computes mass matrix
    mass1d.computeMatrix();

    // prints mass matrix
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < n + 1; j++)
            cout << scientific << mass1d.getMatrixValue(i, j) << " ,  ";
        cout << endl;
    }

    return 0;
}