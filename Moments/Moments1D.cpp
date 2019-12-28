#include <iostream>
#include <cmath>
#include "Moments.h"
#include "JacobiGaussNodes.h"

using std::vector;

BMoment1D::BMoment1D(uint q, uint n, const Element<Element_t::LinearEl> &element)
    : BMoment(q, n, element)
{
    lenMoments = n + 1;
    lenCval = q;
    Bmoment.resize(lenMoments);
    Cval.resize(lenCval);
    intPoints.resize(q);
    intWeights.resize(q);
    assignQuadra();
}

BMoment1D::BMoment1D(const BMoment1D &cp)
    : BMoment(cp.q, cp.n, cp.element)
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
    for (uint k = 0; k < q; k++)
    {
        intPoints[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        intWeights[k] = legendre_w(q, k) * 0.5;
    }
}

void BMoment1D::loadFunctionDef()
{
    for (uint i = 0; i < q; i++)
    {
        Cval[i] = f(intPoints[i]);
    }

    fValSet = true;
}

// returns the vector with the integration points over the object's element, following the moments organization
TPZFMatrix<REAL> BMoment1D::getIntegrationPoints()
{
	TPZFMatrix<REAL> points(q);
	for (uint i = 0; i < q; i++)
	{
		points(i) = intPoints[i];
	}
	return points;
}

// compute the B-moments using the values already assigned in the object
TPZVec<REAL> &BMoment1D::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (!functVal && !fDefSet)
    {
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    }
    else
    {
		std::fill(Bmoment.begin(), Bmoment.end(), 0.0); // zeroes the Bmoment vector
        for (uint i = 0; i < q; i++)
        {
            double xi = intPoints[i];
			double omega = intWeights[i];
            double s = 1 - xi;
            double r = xi / s;
            double w = omega * pow(s, n);
            for (uint alpha = 0; alpha <= n; alpha++)
            {
                // here w equals the Bernstein polynomial of order n,
                // with index alpha, evaluated at the i-th integration node
                // times the i-th integration weight.
                uint p = element.position({alpha});
                Bmoment[p] += w * Cval[i];
                w *= r * ((n - alpha) / (1. + alpha)); // treats the recurrence relation
            }
        }
    }
    return Bmoment;
}
