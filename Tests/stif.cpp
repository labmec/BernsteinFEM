#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

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
        double *Fval = new double[q];
        for (int i = 0; i < q; i++)
            Fval[i] = 1.0; //function( (legendre_xi(q, i) + 1.0) * 0.5 );

        stif.setFunction(Fval);
    }

    // computes stiffness matrix
    stif.compute_matrix();

    // prints stiffness matrix
    for (int i = 0; i < n+1; i++)
    {
        for (int j = 0; j < n+1; j++)
            cout << scientific << stif.getMatrixValue(i, j) << " ,  ";
        cout << endl;
    }
}