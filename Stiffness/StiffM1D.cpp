#include "StiffM.h"

BStiff1D::BStiff1D(int q, int n, const Element<Element_t::LinearEl> &el)
    : BStiff(q, n), BMoment1D(q, 2 * (n - 1), el)
{
    lenStiff = n + 1;
    Matrix.zeros(lenStiff, lenStiff);
}

BStiff1D::~BStiff1D()
{
}

void BStiff1D::computeMatrix()
{
    computeMoments();
    computeBinomials();

    int n = BStiff::n;
    int a = element.getVertices()(0), b = element.getVertices()(1);

    double Const = 1. / BinomialMat.at(n - 1, n - 1) / pow(b - a, n);
    
    for (uint i = 0; i < lenStiff - 1; i++)
    {
        for (uint j = 0; j < lenStiff - 1; j++)
        {
            double w = BinomialMat.at(i, j) * Const;
            w *= BinomialMat.at(n - i - 1, n - j - 1);

            for (uint k = 1; k <= 2; k++)
            {
                for (uint l = 1; l <= 2; l++)
                {
                    uint I = element.position({i + 2 - k}, n);
                    uint J = element.position({j + 2 - l}, n);
                    Matrix.at(I, J) += (n * n) * w * grad(k, l) * Bmoment(i + j);
                }
            }
        }
    }
}