#include "MassM.h"

#define LEN(n) ((n + 1) * (n + 1))

BMass2DQuad::BMass2DQuad(uint q, uint n, Element<Element_t::QuadrilateralEl> const &el)
    : BMass(q, n), BMoment2DQuad(q, 2 * n, el), perm(n, element.getIndexVector())
{
    lenMass = LEN(n);
    Matrix.set_size(lenMass, lenMass);
}

BMass2DQuad::BMass2DQuad(const BMass2DQuad &cp)
    : BMass(cp.BMass::q, cp.BMass::n),
    BMoment2DQuad(cp.BMoment2DQuad::q, cp.BMoment2DQuad::n, cp.element, cp.nb_Array),
    perm(cp.BMass::n, element.getIndexVector())
{
    lenMass = cp.lenMass;
    Matrix = cp.Matrix;
}

BMass2DQuad &BMass2DQuad::operator=(const BMass2DQuad &cp)
{
    if (this != &cp)
    {
        BMass::q = cp.BMass::q;
        BMass::n = cp.BMass::n;
        Matrix = cp.Matrix;
        BinomialMat = cp.BinomialMat;
        lenMass = cp.lenMass;
        lenBinomialMat = cp.lenBinomialMat;
        BMoment2DQuad::q = cp.BMoment2DQuad::q;
        BMoment2DQuad::n = cp.BMoment2DQuad::n;
        lenCval = cp.lenCval;
        lenMoments = cp.lenMoments;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
    }
    return *this;
}

BMass2DQuad::~BMass2DQuad()
{
}

void BMass2DQuad::computeMatrix()
{
    uint n = BMass::n;
    computeMoments();

    double Const = 1.0 / (BinomialMat.at(n, n) * BinomialMat.at(n, n)); // constant due to integration

    // since it is a simple tensor product, this is just like in the 1D case
    // except for indexing
    for (uint a1 = 0; a1 <= n; a1++)
    {
        for (uint a2 = 0; a2 <= n; a2++)
        {
            uint i = perm.getPermutationVector()[a1 * n + a2];
            for (uint b1 = 0; b1 <= n; b1++)
            {
                for (uint b2 = 0; b2 <= n; b2++)
                {
                    double w = Const * BinomialMat.at(a1, b1) * BinomialMat.at(a2, b2);
                    w *= (BinomialMat.at(n - a1, n - b1) * BinomialMat.at(n - a2, n - b2));
                    uint j = perm.getPermutationVector()[b1 * n + b2];
                    Matrix.at(i, j) = w * getBMoment(a1 + b1, a2 + b2, 0);
                }
            }
        }
    }
}