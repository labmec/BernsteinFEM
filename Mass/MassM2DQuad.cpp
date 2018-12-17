#include "MassM.h"

#define LEN(n) ((n + 1) * (n + 1))

BMass2DQuad::BMass2DQuad(int q, int n, Element<Element_t::QuadrilateralEl> el)
    : BMass(q, n), BMoment2DQuad(q, 2 * n, el)
{
    lenMass = LEN(n);

    Matrix.set_size(lenMass, lenMass);
}

BMass2DQuad::~BMass2DQuad()
{
}

void BMass2DQuad::computeMatrix()
{
    int n = BMass::n;
    computeMoments();

    double Const = 1.0 / (BinomialMat.at(n, n) * BinomialMat.at(n, n)); // constant due to integration

    // since it is a simple tensor product, this is just like in the 1D case
    // except for indexing
    for (int a1 = 0; a1 <= n; a1++)
    {
        for (int a2 = 0; a2 <= n; a2++)
        {
            for (int b1 = 0; b1 <= n; b1++)
            {
                for (int b2 = 0; b2 <= n; b2++)
                {
                    double w = Const * BinomialMat.at(a1, b1) * BinomialMat.at(a2, b2);
                    w *= (BinomialMat.at(n - a1, n - b1) * BinomialMat.at(n - a2, n - b2));
                    int i = position(a1, a2, n);
                    int j = position(b1, b2, n);
                    Matrix.at(i, j) = w * getBMoment(a1 + b1, a2 + b2, 0);
                }
            }
        }
    }
}