#include "StiffM.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n) ((n + 1) * (n + 1))

BStiff2DTri::BStiff2DTri(uint32_t q, uint32_t n, const Element<Element_t::TriangularEl> &el)
    : BStiff(q, n), BMoment2DTri(q, 2 * (n - 1), el), normalMat(3, 2)
{
    lenStiff = (n + 1) * (n + 1);
    Matrix.Resize(lenStiff, lenStiff);
}

BStiff2DTri::BStiff2DTri(const BStiff2DTri& cp) :
	BStiff(cp),
	BMoment2DTri(cp),
	normalMat(cp.normalMat)
{}

BStiff2DTri& BStiff2DTri::operator=(const BStiff2DTri& cp) {
    if (this != &cp) {
        BStiff::operator=(cp);
        BMoment2DTri::operator=(cp);
        normalMat = cp.normalMat;
    }
    return *this;
}

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
    const auto norm0 = sqrt(normalMat(0, 0) * normalMat(0, 0) + normalMat(0, 1) * normalMat(0, 1));
    const auto norm1 = sqrt(normalMat(1, 0) * normalMat(1, 0) + normalMat(1, 1) * normalMat(1, 1));
    const auto norm2 = sqrt(normalMat(2, 0) * normalMat(2, 0) + normalMat(2, 1) * normalMat(2, 1));

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
            const auto w1 = Const * BinomialMat(a1, b1);
            for (uint a2 = 0; a2 < n - a1; a2++)
            {
                for (uint b2 = 0; b2 < n - b1; b2++)
                {
                    auto w2 = w1 * BinomialMat(a2, b2) * BinomialMat(n - a1 - a2 - 1, n - b1 - b2 - 1);
                    
                    const uint i = element.position({a1, b1});
                    const uint j = element.position({a2, b2});
                    const uint I = element.position({a1 + 1, b1 + 1});
                    const uint J = element.position({a2 + 1, b2 + 1});
                    const uint I_a = element.position({a1 + 1, b1});
                    const uint J_a = element.position({a2 + 1, b2});
                    const uint I_b = element.position({a1, b1 + 1});
                    const uint J_b = element.position({a2, b2 + 1});

					w2 *= Bmoment[j]; // TODO: check this value

                    const auto n1_n2 = w2 * N[1];
                    const auto n1_n3 = w2 * N[2];
                    const auto n2_n3 = w2 * N[4];

                    Matrix(I, j)     += w2 * N[0];  // n1 . n1
                    Matrix(I_a, J_b) += n1_n2;      // n1 . n2
                    Matrix(I_b, J_a) += n1_n2;      // n2 . n1
                    Matrix(I_a, j)   += n1_n3;      // n1 . n3
                    Matrix(I_b, j)   += n1_n3;      // n3 . n1
                    Matrix(i, J)     += w2 * N[3];  // n2 . n2
                    Matrix(i, J_a)   += n2_n3;      // n2 . n3
                    Matrix(i, J_b)   += n2_n3;      // n3 . n2
                    Matrix(i, j)     += w2 * N[5];  // n3 . n3
                }
            }
        }
    }
}