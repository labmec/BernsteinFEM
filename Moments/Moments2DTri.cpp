#include <iostream>
#include <cmath>
#include <cstdlib>

#include "Moments.h"
#include "JacobiGaussNodes.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n, q) (MAX(n + 1, q) * MAX(n + 1, q))

using std::vector;

// helps indexing quadrature points vectors
uint position_q(uint i, uint j, uint q) { return i * q + j; }
// jacobian determinant of Duffy transform is a multiple of the area of the triangle
double Area2d(double v1[2], double v2[2], double v3[2]);
double Area2d(TPZFMatrix<REAL> &vertices);

BMoment2DTri::BMoment2DTri(uint q, uint n, const Element<Element_t::TriangularEl> &element)
    : BMoment(q, n, element),
      BMomentInter((MAX(n, q - 1) + 1) * (MAX(n, q - 1) + 1))
{
    lenMoments = (n + 1) * (n + 2) * 0.5;
    lenCval = q * q;
    Bmoment.resize(lenMoments);
	Cval.resize(lenCval);
	intPoints.resize(2 * q);
	intWeights.resize(2 * q);
    assignQuadra();
}

BMoment2DTri::BMoment2DTri(const BMoment2DTri &cp)
    : BMoment(cp.q, cp.n, cp.element)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
    intPoints = cp.intPoints;
	intWeights = cp.intWeights;
}

BMoment2DTri &BMoment2DTri::operator=(const BMoment2DTri &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        lenMoments = cp.lenMoments;
        lenCval = cp.lenCval;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
        intPoints = cp.intPoints;
		intWeights = cp.intWeights;
    }
    return *this;
}

BMoment2DTri::~BMoment2DTri() {}

inline
void BMoment2DTri::assignQuadra()
{
    for (uint k = 0; k < q; k++)
    {
        intPoints[k] = (1.0 + jacobi_xi(q, k)) * 0.5;
        intWeights[k] = jacobi_w(q, k) * 0.25;
    }

    for (uint k = q; k < 2*q; k++)
    {
        intPoints[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        intWeights[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment2DTri::loadFunctionDef()
{
    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            uint index_ij = position_q(i, j, q);
            Cval[index_ij] = f(intPoints[i], intPoints[j + q]);
        }
    }
    fValSet = true;
}

// returns the vectors with the integration points (x, y) over the object's element, following the moments organization
// usage: points = getIntegrationPoints(); then
// points(i, 0) == x i-th coordinate
// points(i, 1) == y i-th coordinate
TPZFMatrix<REAL> BMoment2DTri::getIntegrationPoints()
{
    TPZFMatrix<REAL> points(q * q, 2); // vector with q * q * nb_Array elements

    for (uint i = 0; i < q; i++)
    {
		double xi = intPoints[i];
        for (uint j = 0; j < q; j++)
        {
            TPZFMatrix<REAL> coord = {xi, intPoints[j]};
            auto v = element.mapToElement(coord);
            uint index_ij = position_q(i, j, q);
            points(index_ij, 0) = v[0];
            points(index_ij, 1) = v[1];
        }
    }

    return points;
}

//compute the b-moments
TPZVec<REAL> &BMoment2DTri::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'computeMoments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'computeMoments()\'\n";
    else
    {
        uint m = MAX(n + 1, q); // m will be used for indexing
		std::fill(Bmoment.begin(), Bmoment.end(), 0.0); // Zeroes the moments vector

        // compute the function definition into the function values vector
        if (!functVal)
            loadFunctionDef();

        const double scalingConst = Area2d(element.getVertices()); // NEED this scaling constant, as result depends on Area2d(v1,v2,v3)
        
        // convert first index (l=2)
        for (uint i = 0; i < q; i++)
        {
			double xi = intPoints[i];
			double wgt = intWeights[i];

            double s = 1 - xi;
            double r = xi / s;

            double B = wgt * pow(s, n);
            
            for (uint a1 = 0; a1 <= n; a1++)
            {
                for (uint j = 0; j < q; j++)
                {
                    uint index_a1j = position_q(a1, j, m);

                    uint index_ij = position_q(i, j, q);

					BMomentInter[index_a1j] += scalingConst * B * Cval[index_ij];
                }
                B = B * r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=1)
        for (uint i = 0; i < q; i++)
        {
			double xi = intPoints[i];
			double wgt = intWeights[i];

            double s = 1 - xi;
            double r = xi / s;

            for (uint a1 = 0; a1 <= n; a1++)
            {
                double B = wgt * pow(s, n - a1);
                for (uint a2 = 0; a2 <= n - a1; a2++)
                {
                    uint index_a1a2 = element.position({a1, a2});
                    uint index_a1i = position_q(a1, i, m);

					Bmoment[index_a1a2] += B * BMomentInter[index_a1i];
                    B = B * r * (n - a1 - a2) / (1 + a2);
                }
            }
        }
    }
    return Bmoment;
}

// computes area of triangle < v1,v2,v3 >
double Area2d(double v1[2], double v2[2], double v3[2])
{
    double x1 = v1[0];
    double y1 = v1[1];
    double x2 = v2[0];
    double y2 = v2[1];
    double x3 = v3[0];
    double y3 = v3[1];

    return fabs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2.;
}

// computes area of triangle defined by vertices
double Area2d(TPZFMatrix<REAL> &vertices)
{
    double x1 = vertices(0, 0);
    double y1 = vertices(0, 1);
    double x2 = vertices(1, 0);
    double y2 = vertices(1, 1);
    double x3 = vertices(2, 0);
    double y3 = vertices(2, 1);

    return fabs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2.;
}