#include <cmath>

#include "Elem.h"
#include "JacobiGaussNodes.h"

BElement2DQuad::BElement2DQuad(int q, int n)
    : MassMat(q, n), StiffMat(q, n), LoadVec(q, n),
      ElMat((n + 1) * (n + 1), (n + 1) * (n + 1), arma::fill::none),
      BBVector(n + 1, arma::fill::none),
      QuadVector(q, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    length = (n + 1) * (n + 1);
}

BElement2DQuad::~BElement2DQuad()
{
}

void BElement2DQuad::setQuad(double v1[2], double v2[2], double v3[2], double v4[2])
{
    MassMat.setQuadrilateral(v1, v2, v3, v4);
    //ConvecMat.setQuadrilateral(v1, v2, v3, v4);
    StiffMat.setQuadrilateral(v1, v2, v3, v4);
    LoadVec.setQuadrilateral(v1, v2, v3, v4);
}

arma::vec BElement2DQuad::evaluate()
{
    // initialize QuadVector with 0's
    for (int i = 0; i < q * q; i++)
        QuadVector(i) = 0.0;

    // EvalStep routine, same as the 1D with tensor product
    for (int i = 0; i < q; i++)
    {
        double xi1 = legendre_xi(q, i);
        double s1 = 1 - xi1;
        double r1 = xi1 / s1;
        double w1 = pow(s1, n);
        for (int j = 0; j < q; j++)
        {
            double xi2 = legendre_xi(q, j);
            double s2 = 1 - xi2;
            double r2 = xi2 / s2;
            double w2 = pow(s2, n);

            for (int a = 0; a < n + 1; a++)
            {
                for (int b = 0; b < n + 1; b++)
                {
                    QuadVector[n * i + j] += w2 * BBVector[n * i + j];
                    w2 *= r2 * (n / (1. + b));
                }
                w1 *= r1 * (n / 1. + a);
            }
        }
    }
    return QuadVector;
}

void BElement2DQuad::makeSystem()
{
    MassMat.compute_matrix();
    //ConvecMat.compute_matrix();
    StiffMat.compute_matrix();
    LoadVec.compute_moments();

    for (int i = 0; i < length; i++)
        for (int j = 0; j < length; j++)
            ElMat(i, j) = MassMat.getMatrixValue(i, j) + StiffMat.getMatrixValue(i, j);
}