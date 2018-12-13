#include <iostream>
#include <cmath>
#include <cstdlib>
using std::cin;
using std::cout;
using std::endl;

#include "Moments.h"
#include "JacobiGaussNodes.h"

inline
void BMoment1D::assignQuadra()
{
    double *x = quadraWN.colptr(1);
    double *w = quadraWN.colptr(0);

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment1D::loadFunctionDef()
{
    arma::mat points = getIntegrationPoints();
    for (int i = 0; i < q; i++)
    {
        for (int el = 0; el < nb_Array; el++)
            Cval.at(i, el) = f(points(i));
    }

    fValSet = true;
}

// returns the vector with the integration points over the object's element, following the moments organization
inline
arma::mat BMoment1D::getIntegrationPoints()
{
    arma::mat points(q, nb_Array, arma::fill::none); // vector with q elements
    double a = element.getVertices()(0);
    double b = element.getVertices()(1);

    for(int i = 0; i < q; i++)
    {
        points.at(i) = (legendre_xi(q, i) + (a + 1)) * ((b - a) / 2.);
    }

    return points;
}

// compute the B-moments using the values already assigned in the object
inline
void BMoment1D::computeMoments()
{
    if (functVal == 0 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 1 && !fValSet)
    {
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
        loadFunctionDef();
    }
    else
    {
        Bmoment.zeros();
        for (int i = 0; i < q; i++)
        {
            double xi = quadraWN.at(i, 1);
            double omega = quadraWN.at(i, 0);
            double s = 1 - xi;
            double r = xi / s;
            double w = omega * pow(s, n);
            for (int alpha = 0; alpha <= n; alpha++)
            {
                // here w equals the Bernstein polynomial of order n,
                // with index alpha, evaluated at the i-th integration node
                // times the i-th integration weight.
                int p = position(alpha, n);
                for (int el = 0; el < nb_Array; el++)
                    Bmoment.at(p, el) += w * Cval.at(i, el);
                w *= r * ((n - alpha) / (1. + alpha)); // treats the recurrence relation
            }
        }
    }
}
