#include <cmath>

#include "Elem.h"
#include "JacobiGaussNodes.h"

BElement2DTri::BElement2DTri(int q, int n)
    : MassMat(q, n), StiffMat(q, n), LoadVec(q, n),
      ElMat((n + 1) * (n + 1), (n + 1) * (n + 1), arma::fill::none),
      BBVector((n + 1) * (n + 1), arma::fill::none),
      QuadVector(q * q, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    //length = ((n + 1) * (n + 2)) / 2;
    length = (n + 1) * (n + 1); // accounts for positioning
}

BElement2DTri::~BElement2DTri()
{
}

void BElement2DTri::setTriangle(double v1[2], double v2[2], double v3[2])
{
    MassMat.setTriangle(v1, v2, v3);
    //ConvecMat.setTriangle(v1, v2, v3);
    StiffMat.setTriangle(v1, v2, v3);
    LoadVec.setTriangle(v1, v2, v3);
}

arma::vec BElement2DTri::evaluate()
{
    int Length = ((n + 1) * (n + 2)) / 2;
    arma::vec QuadInter(q * q, arma::fill::zeros);

    // convert first index
    for (int i = 0; i < q; i++)
    {
        double xi = legendre_xi(q, i);
        double s = 1 - xi;
        double r = xi / s;

        for (int a1 = 0; a1 < Length; a1++)
        {
            double w = pow(s, n - a1);
            for (int a2 = 0; a2 < n - a1; a2++)
            {
                QuadInter[i] += w * BBVector(BMoment2DTri::position(a1, a2, n));
                w *= r * (n - a1 - a2) / (1. + a2);
            }
        }
    }

    // convert second index
    for (int i = 0; i < q; i++)
    {
        double xi = jacobi_xi(q, i);
        double s = 1 - xi;
        double r = xi / s;
        double w = pow(s, n);
        for (int a1 = 0; a1 < Length; a1++)
        {
            for (int i2 = 0; i2 < q; i2++)
            {
                QuadVector[i] += w * QuadInter[i2];
            }
            w *= r * (n - a1 / (1. + a1));
        }
    }

    return QuadVector;
}

void BElement2DTri::makeSystem()
{
    MassMat.compute_matrix();
    //ConvecMat.compute_matrix();
    StiffMat.compute_matrix();
    LoadVec.compute_moments();

    for (int i = 0; i < length; i++)
        for (int j = 0; j < length; j++)
            ElMat(i, j) = MassMat.getMatrixValue(i, j) + StiffMat.getMatrixValue(i, j);
}