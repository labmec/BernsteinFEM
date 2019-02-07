#include "StiffM.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n) ((n + 1) * (n + 1))

using namespace arma;

// functions implemented in another file (Moments/Moments2DTri.cpp)
extern double Area2d(double v1[2], double v2[2], double v3[2]);
extern double Area2d(const arma::mat &vertices);

BStiff2DTri::BStiff2DTri(uint q, uint n, const Element<Element_t::TriangularEl> &el)
    : BStiff(q, n), BMoment2DTri(q, 2 * (n - 1), el), normalMat(3, 2, arma::fill::none)
{
    lenStiff = (n + 1) * (n + 1);
    Matrix.set_size(lenStiff, lenStiff);
}

BStiff2DTri::~BStiff2DTri() { }

void BStiff2DTri::compute_normals()
{
    auto vertices = element.getVertices();
    normalMat(0, 0) = vertices(1, 1) - vertices(2, 1);
    normalMat(0, 1) = vertices(2, 0) - vertices(1, 0);
    normalMat(1, 0) = vertices(2, 1) - vertices(0, 1);
    normalMat(1, 1) = vertices(0, 0) - vertices(2, 0);
    normalMat(2, 0) = vertices(0, 1) - vertices(1, 1);
    normalMat(2, 1) = vertices(1, 0) - vertices(0, 0);

    // normalize_normals();
    // the normals computed this way have the same length as the opposite side of the vertex
    // which is a property we're gonna make use of, so there's no need to normalize
}

void BStiff2DTri::normalize_normals()
{
    double norm0 = sqrt(normalMat.at(0, 0) * normalMat.at(0, 0) + normalMat.at(0, 1) * normalMat.at(0, 1));
    double norm1 = sqrt(normalMat.at(1, 0) * normalMat.at(1, 0) + normalMat.at(1, 1) * normalMat.at(1, 1));
    double norm2 = sqrt(normalMat.at(2, 0) * normalMat.at(2, 0) + normalMat.at(2, 1) * normalMat.at(2, 1));

    normalMat.at(0, 0) = normalMat.at(0, 0) / norm0;
    normalMat.at(0, 1) = normalMat.at(0, 1) / norm0;
    normalMat.at(1, 0) = normalMat.at(1, 0) / norm1;
    normalMat.at(1, 1) = normalMat.at(1, 1) / norm1;
    normalMat.at(2, 0) = normalMat.at(2, 0) / norm2;
    normalMat.at(2, 1) = normalMat.at(2, 1) / norm2;
}

void BStiff2DTri::compute_normals_products(vec &N)
{
    if (N.size() < 6) {
        N.resize(6);
    }
    N[0] = normalMat.at(0, 0) * normalMat.at(0, 0) + normalMat.at(0, 1) * normalMat.at(0, 1); // n1 . n1
    N[1] = normalMat.at(0, 0) * normalMat.at(1, 0) + normalMat.at(0, 1) * normalMat.at(1, 1); // n1 . n2
    N[2] = normalMat.at(0, 0) * normalMat.at(2, 0) + normalMat.at(0, 1) * normalMat.at(2, 1); // n1 . n3
    N[3] = normalMat.at(1, 0) * normalMat.at(1, 0) + normalMat.at(1, 1) * normalMat.at(1, 1); // n2 . n2
    N[4] = normalMat.at(1, 0) * normalMat.at(2, 0) + normalMat.at(1, 1) * normalMat.at(2, 1); // n2 . n3
    N[5] = normalMat.at(2, 0) * normalMat.at(2, 0) + normalMat.at(2, 1) * normalMat.at(2, 1); // n3 . n3
}

void BStiff2DTri::computeMatrix()
{
    vec N(6, fill::none);
    computeMoments();           // computes the moments necessary for the computation of the stiffness matrix
    compute_normals();           // computes the normal vectors, used in the gradient
    compute_normals_products(N); // computes the product of each combination of the normals
    // transform_BmomentC_Stiff2d(Moments, normalMat);

    Matrix.zeros();
    uint n = BStiff::n;
    double area = Area2d(element.getVertices());
    double Const = n * n / 4. / area / area / BinomialMat(n - 1, n - 1); // taking account of scaling between normals and gradients

    
    // com
    for (uint a1 = 0; a1 < n; a1++)
    {
        for (uint b1 = 0; b1 < n; b1++)
        {
            double w1 = Const * BinomialMat.at(a1, b1);
            for (uint a2 = 0; a2 < n - a1; a2++)
            {
                for (uint b2 = 0; b2 < n - b1; b2++)
                {
                    double w2 = w1 * BinomialMat.at(a2, b2) * BinomialMat.at(n - a1 - a2 - 1, n - b1 - b2 - 1);
                    w2 *= getBMoment(a1 + b1, a2 + b2, 0);

                    uint i = element.position({a1, b1}, n);
                    uint j = element.position({a2, b2}, n);
                    uint I = element.position({a1 + 1, b1 + 1}, n);
                    uint J = element.position({a2 + 1, b2 + 1}, n);
                    uint I_a = element.position({a1 + 1, b1}, n);
                    uint J_a = element.position({a2 + 1, b2}, n);
                    uint I_b = element.position({a1, b1 + 1}, n);
                    uint J_b = element.position({a2, b2 + 1}, n);

                    double n1n2 = w2 * N[1];
                    double n1n3 = w2 * N[2];
                    double n2n3 = w2 * N[4];

                    Matrix.at(I, j)     += w2 * N[0]; // n1 . n1
                    Matrix.at(I_a, J_b) += n1n2;      // n1 . n2
                    Matrix.at(I_b, J_a) += n1n2;      // n2 . n1
                    Matrix.at(I_a, j)   += n1n3;      // n1 . n3
                    Matrix.at(I_b, j)   += n1n3;      // n3 . n1
                    Matrix.at(i, J)     += w2 * N[3]; // n2 . n2
                    Matrix.at(i, J_a)   += n2n3;      // n2 . n3
                    Matrix.at(i, J_b)   += n2n3;      // n3 . n2
                    Matrix.at(i, j)     += w2 * N[5]; // n3 . n3
                }
            }
        }
    }
}
