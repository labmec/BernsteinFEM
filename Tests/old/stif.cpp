/** stif.cpp
 * 
 * This file makes a test for the computation
 * of the finite element Stiffness Matrix for
 * 1D elements using Bernstein polynomials
 */

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <armadillo>

#include "StiffM.h"
#include "JacobiGaussNodes.h"

using namespace std;

double function(double x)
{
    return sin(M_PI * x);
}

int main()
{
    int n = 4;      // basis order
    int q = 2 * n;  // number of integration points

    // geometry element (default: [0,1])
    // Element<Element_t::LinearEl> el({0, 1});

    // Mass matrix object
    // BStiff1D stif(q, n, el);
    BStiff1D stif(q, n);

    // sets function values (1)
    {
        arma::vec Fval(q, arma::fill::ones);

        /* for (int i = 0; i < q; i++)
            Fval[i] = 1.0; //function( (legendre_xi(q, i) + 1.0) * 0.5 ); */

        stif.setFunctionValues(Fval);
    }

    // computes stiffness matrix
    stif.computeMatrix();

    // prints stiffness matrix
    cout << "{";
    for (int i = 0; i < (n+1); i++)
    {
        cout << "{";
        for (int j = 0; j < (n+1); j++)
            if (j < n)
                cout << stif.getMatrixValue(i, j) << ", ";
            else   
                cout << stif.getMatrixValue(i, j);
        if (i < n)
            cout << "}," << endl;
        else
            cout << "}";
    }
    cout << "}" << endl;

	return 0;
}