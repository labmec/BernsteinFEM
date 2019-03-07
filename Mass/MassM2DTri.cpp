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
    : BMass(q, n), BMoment2DTri(q, 2 * n, el)
{
    lenMass = ((n + 1) * (n + 1));
    Matrix.set_size(lenMass, lenMass);
}

BMass2DTri::BMass2DTri(const BMass2DTri &cp)
    : BMass(cp.BMass::q, cp.BMass::n), BMoment2DTri(cp.BMoment2DTri::q, cp.BMoment2DTri::n, cp.element)
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

// this method was based on the code by M. Ainsworth
// G. Andriamro and O. Davydov
// void BMass2DTri::compute_matrix()
// {
//     compute_moments();
//     compute_binomials();

//     double Const = 1.0 / BinomialMat.at(n, n);

//     int M = MAX(2 * n, q - 1);

//     int iMu = 0;
//     for (int mu0 = n; mu0 >= 0; mu0--)
//     {
//         int iEta = 0;
//         for (int eta0 = n; eta0 >= 0; eta0--)
//         {
//             int iMuEta0 = BMoment2DTri::position(mu0, eta0, M); // the Bmoments are indexed w.r.t position2d2(. , MAX(2n,q-1) )
//             double w1 = Const * BinomialMat.at(mu0, eta0);

//             int mu2 = 0;
//             for (int mu1 = n - mu0; mu1 >= 0; mu1--, mu2++, iMu++)
//             {
//                 int eta2 = 0;
//                 for (int eta1 = n - eta0; eta1 >= 0; eta1--, eta2++, iEta++)
//                 {
//                     double w2 = w1 * BinomialMat.at(mu1, eta1) * BinomialMat.at(mu2, eta2);

//                     int iMuEta = iMuEta0 + mu1 + eta1;

//                     Matrix.at(iMu, iEta) = w2 * get_bmoment(iMuEta);
//                 }
//                 if (mu1 > 0)
//                     iEta -= eta2;
//             }
//             if (eta0 > 0)
//                 iMu -= mu2;
//         }
//     }
// }

void BMass2DTri::computeMatrix()
{
    uint n = BMass::n;
    computeMoments();

    double Const = 1.0 / BinomialMat(n, n);

    for (uint a1 = 0; a1 <= n; a1++) {
        for (uint b1 = 0; b1 <= n; b1++) {
            double w1 = Const * BinomialMat(a1, b1);
            for (uint a2 = 0; a2 <= n - a1; a2++){
                for (uint b2 = 0; b2 <= n - b1; b2++) {
                    double w2 = w1 * BinomialMat(a2, b2);
                    uint i = element.position({a1, a2});
                    uint j = element.position({b1, b2});
                    Matrix(i, j) = w2 * getBMoment(a1 + b1, a2 + b2, 0);
                } 
            }
        }
    }
}