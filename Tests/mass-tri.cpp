#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>

#include <armadillo>

#include "MassM.h"
#include "JacobiGaussNodes.h"

using namespace std;

double f(double x, double y) {
    return 2 - sin(x*y);
}

int main()
{
    int n = 4;
    int q = 2 * n;

    BMass2DTri mass_tri(q, n);

    { // sets Triangle
        double v1[2] = {0, 0};
        double v2[2] = {1, 0};
        double v3[2] = {0, 1};
        mass_tri.setTriangle(v1, v2, v3);
    }

    // sets function values (1)
    {
        arma::vec Fval(q * q, arma::fill::ones);
        mass_tri.setFunction(Fval);
    }

    // computes mass matrix
    mass_tri.compute_matrix();

    // for (int a1 = 0; a1 < n; a1++) {
    //     for (int a2 = 0; a2 < n - a1; a2 ++)
    //         cout << "Bmoment[" << a1 << ", " << a2 << ", " << n - a1 - a2 << "] : " << mass_tri.get_bmoment(a1, a2, 2 * n);
    // }

    // prints mass matrix
    int iMu = 0;
    for (int mu0 = n; mu0 >= 0; mu0--)
    {
        int mu2 = 0;
        for (int mu1 = n - mu0; mu1 >= 0; mu1--, mu2++, iMu++)
        {
            int iEta = 0;
            for (int eta0 = n; eta0 >= 0; eta0--)
            {
                int eta2 = 0;
                for (int eta1 = n - eta0; eta1 >= 0; eta1--, eta2++, iEta++)
                {
                    cout << scientific << mass_tri.getMatrixValue(iMu, iEta) << ", ";
                }
            }
            cout << endl;
        }
    }

    return 0;
}