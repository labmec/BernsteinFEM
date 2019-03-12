// TODO: test
#include "MassM.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n) ((n + 1) * (n + 1)) // (((n+1) * (n+2)) / 2)
#ifdef LENB
#undef LENB
#endif
#define LENB(n) (2 * n + 2)

BMass2DTri::BMass2DTri(uint q, uint n, const Element<Element_t::TriangularEl> &el)
    : BMass(q, n), BMoment2DTri(q, 2 * n, el), perm(PermutationPool<Element_t::TriangularEl>::GetPermutation(n))
{
    lenMass = ((n + 1) * (n + 2) / 2);
    Matrix.set_size(lenMass, lenMass);
}

BMass2DTri::BMass2DTri(const BMass2DTri &cp)
    : BMass(cp.BMass::q, cp.BMass::n), BMoment2DTri(cp.BMoment2DTri::q, cp.BMoment2DTri::n, cp.element), perm(PermutationPool<Element_t::TriangularEl>::GetPermutation(cp.BMass::n))
{
    lenMass = cp.lenMass;
    Matrix = cp.Matrix;
}

BMass2DTri &BMass2DTri::operator=(const BMass2DTri &cp)
{
    if (this != &cp)
    {
        BMass::q = cp.BMass::q;
        BMass::n = cp.BMass::n;
        Matrix = cp.Matrix;
        BinomialMat = cp.BinomialMat;
        lenMass = cp.lenMass;
        lenBinomialMat = cp.lenBinomialMat;
        BMoment2DTri::q = cp.BMoment2DTri::q;
        BMoment2DTri::n = cp.BMoment2DTri::n;
        lenCval = cp.lenCval;
        lenMoments = cp.lenMoments;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
    }
    return *this;
}

BMass2DTri::~BMass2DTri()
{
}

void BMass2DTri::computeMatrix()
{
    uint n = BMass::n;
    computeMoments();

    double Const = 1.0 / BinomialMat(n, n);

    for (uint a1 = 0; a1 <= n; a1++) {
        for (uint b1 = 0; b1 <= n; b1++) {
            double w1 = Const * BinomialMat(a1, b1);
            for (uint a2 = 0; a2 <= n; a2++) {
                uint i = perm.getPermutationVector()[a1 * n + a2];
                for (uint b2 = 0; b2 <= n; b2++) {
                    double w2 = w1 * BinomialMat(a2, b2);
                    uint j = perm.getPermutationVector()[b1 * n + b2];
                    Matrix(i, j) = w2 * getBMoment(a1 + b1, a2 + b2, 0);
                } 
            }
        }
    }
}