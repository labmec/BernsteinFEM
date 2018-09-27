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
    int n = 4;
    int q = 2 * n;

    BStiff1D stif(q, n);

    // sets interval
    {
        double a = 0.0;
        double b = 1.0;
        stif.setInterval(a, b);
    }
    

    // sets function values (1)
    {
        arma::vec Fval(q, arma::fill::ones);

        /* for (int i = 0; i < q; i++)
            Fval[i] = 1.0; //function( (legendre_xi(q, i) + 1.0) * 0.5 ); */

        stif.setFunction(Fval);
    }

    // computes stiffness matrix
    stif.compute_matrix();

    // prints stiffness matrix
    cout << "{";
    for (int i = 0; i < (n+1) * (n+1); i++)
    {
        cout << "{";
        for (int j = 0; j < (n+1) * (n+1); j++)
            if (j < (n+1)*(n+1) - 1)
                cout << stif.getMatrixValue(i, j) << ", ";
            else   
                cout << stif.getMatrixValue(i, j);
        if (i < (n+1)*(n+1) - 1)
            cout << "},";
        else
            cout << "}";
    }
    cout << "}" << endl;

	return 0;
}