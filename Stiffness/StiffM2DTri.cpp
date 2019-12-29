#include "StiffM.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n) ((n + 1) * (n + 1))

BStiff2DTri::BStiff2DTri(uint q, uint n, const Element<Element_t::TriangularEl> &el)
    : BStiff(q, n), BMoment2DTri(q, 2 * (n - 1), el), normalMat(3, 2)
{
    lenStiff = (n + 1) * (n + 1);
    Matrix.Resize(lenStiff, lenStiff);
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
    double norm0 = sqrt(normalMat(0, 0) * normalMat(0, 0) + normalMat(0, 1) * normalMat(0, 1));
    double norm1 = sqrt(normalMat(1, 0) * normalMat(1, 0) + normalMat(1, 1) * normalMat(1, 1));
    double norm2 = sqrt(normalMat(2, 0) * normalMat(2, 0) + normalMat(2, 1) * normalMat(2, 1));

    normalMat(0, 0) = normalMat(0, 0) / norm0;
    normalMat(0, 1) = normalMat(0, 1) / norm0;
    normalMat(1, 0) = normalMat(1, 0) / norm1;
    normalMat(1, 1) = normalMat(1, 1) / norm1;
    normalMat(2, 0) = normalMat(2, 0) / norm2;
    normalMat(2, 1) = normalMat(2, 1) / norm2;
}

void BStiff2DTri::compute_normals_products(TPZVec<REAL> &N)
{
    if (N.size() < 6) {
        N.resize(6);
    }
    N[0] = normalMat(0, 0) * normalMat(0, 0) + normalMat(0, 1) * normalMat(0, 1); // n1 . n1
    N[1] = normalMat(0, 0) * normalMat(1, 0) + normalMat(0, 1) * normalMat(1, 1); // n1 . n2
    N[2] = normalMat(0, 0) * normalMat(2, 0) + normalMat(0, 1) * normalMat(2, 1); // n1 . n3
    N[3] = normalMat(1, 0) * normalMat(1, 0) + normalMat(1, 1) * normalMat(1, 1); // n2 . n2
    N[4] = normalMat(1, 0) * normalMat(2, 0) + normalMat(1, 1) * normalMat(2, 1); // n2 . n3
    N[5] = normalMat(2, 0) * normalMat(2, 0) + normalMat(2, 1) * normalMat(2, 1); // n3 . n3
}

void BStiff2DTri::computeMatrix()
{
    TPZVec<REAL> N(6);
    computeMoments();            // computes the moments necessary for the computation of the stiffness matrix
    compute_normals();           // computes the normal vectors, used in the gradient
    compute_normals_products(N); // computes the product of each combination of the normals
    // transform_BmomentC_Stiff2d(Moments, normalMat);

    Matrix.Zero();
    const auto n = BStiff::n;
    const auto area = triangle_area(element.getVertices());
    const auto Const = n * n / 4. / area / area / BinomialMat(n - 1, n - 1); // taking account of scaling between normals and gradients
    
    // compute stiffness matrix
    for (uint a1 = 0; a1 < n; a1++)
    {
        for (uint b1 = 0; b1 < n; b1++)
        {
            double w1 = Const * BinomialMat(a1, b1);
            for (uint a2 = 0; a2 < n - a1; a2++)
            {
                for (uint b2 = 0; b2 < n - b1; b2++)
                {
                    double w2 = w1 * BinomialMat(a2, b2) * BinomialMat(n - a1 - a2 - 1, n - b1 - b2 - 1);
                    
                    uint i = element.position({a1, b1});
                    uint j = element.position({a2, b2});
                    uint I = element.position({a1 + 1, b1 + 1});
                    uint J = element.position({a2 + 1, b2 + 1});
                    uint I_a = element.position({a1 + 1, b1});
                    uint J_a = element.position({a2 + 1, b2});
                    uint I_b = element.position({a1, b1 + 1});
                    uint J_b = element.position({a2, b2 + 1});

					w2 *= Bmoment[j]; // TODO: check this value

                    double n1n2 = w2 * N[1];
                    double n1n3 = w2 * N[2];
                    double n2n3 = w2 * N[4];

                    Matrix(I, j)     += w2 * N[0]; // n1 . n1
                    Matrix(I_a, J_b) += n1n2;      // n1 . n2
                    Matrix(I_b, J_a) += n1n2;      // n2 . n1
                    Matrix(I_a, j)   += n1n3;      // n1 . n3
                    Matrix(I_b, j)   += n1n3;      // n3 . n1
                    Matrix(i, J)     += w2 * N[3]; // n2 . n2
                    Matrix(i, J_a)   += n2n3;      // n2 . n3
                    Matrix(i, J_b)   += n2n3;      // n3 . n2
                    Matrix(i, j)     += w2 * N[5]; // n3 . n3
                }
            }
        }
    }
}
