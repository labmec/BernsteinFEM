#include <iostream>
#include <cmath>
#include "Moments.h"
#include "JacobiGaussNodes.h"

BMoment1D::BMoment1D(uint q, uint n, const Element<Element_t::LinearEl> &element, uint nb_Array)
    : BMoment(q, n, element, nb_Array)
{
    lenMoments = n + 1;
    lenCval = q;
    Bmoment.set_size(lenMoments, nb_Array);
    quadraWN.set_size(q);
    Cval.set_size(lenCval, nb_Array);
    assignQuadra();
}

BMoment1D::BMoment1D(const BMoment1D &cp)
    : BMoment(cp.q, cp.n, cp.element, cp.nb_Array)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
}

BMoment1D &BMoment1D::operator=(const BMoment1D &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        nb_Array = cp.nb_Array;
        lenMoments = cp.lenMoments;
        lenCval = cp.lenCval;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
    }
    return *this;
}

BMoment1D::~BMoment1D() {}

inline
void BMoment1D::assignQuadra()
{
    double *x = quadraWN.colptr(1);
    double *w = quadraWN.colptr(0);

    for (uint k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment1D::loadFunctionDef()
{
    arma::mat points = getIntegrationPoints();
    for (uint i = 0; i < q; i++)
    {
        for (uint el = 0; el < nb_Array; el++)
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

    for(uint i = 0; i < q; i++)
    {
        points.at(i) = (legendre_xi(q, i) + (a + 1)) * ((b - a) / 2.);
    }

    return points;
}

// compute the B-moments using the values already assigned in the object
arma::mat &BMoment1D::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (!functVal && !fDefSet)
    {
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    }
    else
    {
        Bmoment.zeros();
        for (uint i = 0; i < q; i++)
        {
            double xi = quadraWN.at(i, 1);
            double omega = quadraWN.at(i, 0);
            double s = 1 - xi;
            double r = xi / s;
            double w = omega * pow(s, n);
            for (uint alpha = 0; alpha <= n; alpha++)
            {
                // here w equals the Bernstein polynomial of order n,
                // with index alpha, evaluated at the i-th integration node
                // times the i-th integration weight.
                uint p = element.position({alpha});
                for (uint el = 0; el < nb_Array; el++)
                    Bmoment.at(p, el) += w * Cval.at(i, el);
                w *= r * ((n - alpha) / (1. + alpha)); // treats the recurrence relation
            }
        }
    }
    return Bmoment;
}
