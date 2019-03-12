#include <iostream>
#include <cmath>
#include <cstdlib>

#include "Moments.h"
#include "JacobiGaussNodes.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n, q) (MAX(n + 1, q) * MAX(n + 1, q))

// helps indexing quadrature points vectors
uint position_q(uint i, uint j, uint q) { return i * q + j; }
// jacobian determinant of Duffy transform is a multiple of the area of the triangle
double Area2d(double v1[2], double v2[2], double v3[2]);
double Area2d(const arma::mat &vertices);

BMoment2DTri::BMoment2DTri(uint q, uint n, const Element<Element_t::TriangularEl> &element, uint nb_Array)
    : BMoment(q, n, element, nb_Array),
      BMomentInter((MAX(n, q - 1) + 1) * (MAX(n, q - 1) + 1), nb_Array, arma::fill::zeros)
{
    lenMoments = (n + 1) * (n + 2) * 0.5;
    lenCval = q * q;
    Bmoment.set_size(lenMoments, nb_Array);
    assignQuadra();
}

BMoment2DTri::BMoment2DTri(const BMoment2DTri &cp)
    : BMoment(cp.q, cp.n, cp.element, cp.nb_Array)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
    quadraWN = cp.quadraWN;
}

BMoment2DTri &BMoment2DTri::operator=(const BMoment2DTri &cp)
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
        quadraWN = cp.quadraWN;
    }
    return *this;
}

BMoment2DTri::~BMoment2DTri() {}

inline
void BMoment2DTri::assignQuadra()
{
    quadraWN.zeros(q, 4);
    double *x = quadraWN.colptr(1); // x is pointer which points to the address l_x, thus the effective MODIFICATION of l_x entries
    double *w = quadraWN.colptr(0);

    for (uint k = 0; k < q; k++)
    {
        x[k] = (1.0 + jacobi_xi(q, k)) * 0.5;
        w[k] = jacobi_w(q, k) * 0.25;
    }

    x = quadraWN.colptr(3);
    w = quadraWN.colptr(2);

    for (uint k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment2DTri::loadFunctionDef()
{
    arma::mat points(getIntegrationPoints());
    Cval.set_size(lenCval, nb_Array);
    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            uint index_ij = position_q(i, j, q);
            Cval(index_ij, 0) = f(points(i, 0), points(j, 1));
        }
    }
    fValSet = true;
}

// returns the vectors with the integration points (x, y) over the object's element, following the moments organization
// usage: points = getIntegrationPoints(); then
// points(i, 0) == x i-th coordinate
// points(i, 1) == y i-th coordinate
inline
arma::mat BMoment2DTri::getIntegrationPoints()
{
    arma::mat points(q * q * nb_Array, 2); // vector with q * q * nb_Array elements

    for (uint i = 0; i < q; i++)
    {
        double xi = quadraWN(i, 1);
        for (uint j = 0; j < q; j++)
        {
            arma::mat coord = {xi, quadraWN(j, 3)};
            arma::mat v = element.mapToElement(coord);
            uint index_ij = position_q(i, j, q);
            points(index_ij, 0) = v[0];
            points(index_ij, 1) = v[1];
        }
    }

    return points;
}

//compute the b-moments
arma::mat &BMoment2DTri::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'computeMoments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'computeMoments()\'\n";
    else
    {
        uint m = MAX(n + 1, q); // m will be used for indexing
        Bmoment.zeros();

        // compute the function definition into the function values vector
        if (!functVal)
            loadFunctionDef();

        const double scalingConst = Area2d(element.getVertices()); // NEED this scaling constant, as result depends on Area2d(v1,v2,v3)
        
        // convert first index (l=2)
        for (uint i = 0; i < q; i++)
        {
            double xi = quadraWN.at(i, 1);
            double wgt = quadraWN.at(i, 0);

            double s = 1 - xi;
            double r = xi / s;

            double B = wgt * pow(s, n);
            
            for (uint a1 = 0; a1 <= n; a1++)
            {
                for (uint j = 0; j < q; j++)
                {
                    uint index_a1j = position_q(a1, j, m);

                    uint index_ij = position_q(i, j, q);

                    for (uint ell = 0; ell < nb_Array; ell++)
                    {
                        BMomentInter.at(index_a1j, ell) += scalingConst * B * Cval.at(index_ij, ell);
                    }
                }
                B = B * r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=1)
        for (uint i = 0; i < q; i++)
        {
            double xi = quadraWN.at(i, 3);
            double wgt = quadraWN.at(i, 2);

            double s = 1 - xi;
            double r = xi / s;

            for (uint a1 = 0; a1 <= n; a1++)
            {
                double B = wgt * pow(s, n - a1);
                for (uint a2 = 0; a2 <= n - a1; a2++)
                {
                    uint index_a1a2 = element.position({a1, a2});
                    uint index_a1i = position_q(a1, i, m);

                    for (uint ell = 0; ell < nb_Array; ell++)
                    {
                        Bmoment(index_a1a2, ell) += B * BMomentInter(index_a1i, ell);
                    }
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

    return abs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2;
}

// computes area of triangle defined by vertices
double Area2d(const arma::mat &vertices)
{
    double x1 = vertices(0, 0);
    double y1 = vertices(0, 1);
    double x2 = vertices(1, 0);
    double y2 = vertices(1, 1);
    double x3 = vertices(2, 0);
    double y3 = vertices(2, 1);

    return fabs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2.;
}