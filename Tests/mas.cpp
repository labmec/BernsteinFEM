#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include "MassM.h"
#include "JacobiGaussNodes.h"

using namespace std;

int main()
{
    int n = 4;
    int q = 2 * n;

    // interval variables
    double a = 0.0;
    double b = 1.0;

    BMass1D mass1d(q, n);

    mass1d.setInterval(a, b);

    // sets function values (1)
    {
        double *Fval = new double[q];
        for (int i = 0; i < q; i++)
            Fval[i] = 1.0;

        mass1d.setFunction(Fval);
    }

    // computes mass matrix
    mass1d.compute_matrix();

    // prints mass matrix
    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < n + 1; j++)
            cout << scientific << mass1d.getMatrixValue(i, j) << " ,  ";
        cout << endl;
    }

    return 0;
}