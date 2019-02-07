#include <iostream>
using std::cin;
using std::cout;
using std::endl;

#include "MassM.h"

// maybe trade q to 2*q in the base class constructor
BMass1D::BMass1D(uint q, uint n, const Element<Element_t::LinearEl> &el)
    : BMass(q, n), BMoment1D(q, 2 * n, el)
{
    lenMass = n + 1;
    Matrix.set_size(lenMass, lenMass);
}

BMass1D::BMass1D(const BMass1D &cp)
    : BMass(cp.BMass::q, cp.BMass::n), BMoment1D(cp.BMoment1D::q, cp.BMoment1D::n, cp.element, cp.nb_Array)
{
    lenMass = cp.lenMass;
    Matrix = cp.Matrix;
}

BMass1D &BMass1D::operator=(const BMass1D &cp)
{
    if (this != &cp)
    {
        BMass::q = cp.BMass::q;
        BMass::n = cp.BMass::n;
        Matrix = cp.Matrix;
        BinomialMat = cp.BinomialMat;
        lenMass = cp.lenMass;
        lenBinomialMat = cp.lenBinomialMat;
        BMoment1D::q = cp.BMoment1D::q;
        BMoment1D::n = cp.BMoment1D::n;
        lenCval = cp.lenCval;
        lenMoments = cp.lenMoments;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
    }
    return *this;
}

BMass1D::~BMass1D()
{
}

void BMass1D::computeMatrix()
{
    uint n = BMass::n;
    computeMoments();

    double Const = 1.0 / BinomialMat.at(n, n);

    for (uint i = 0; i < lenMass; i++)
    {
        uint I = element.position({i}, n);
        for (uint j = 0; j < lenMass; j++)
        {
            double binom = Const * BinomialMat.at(i, j) * BinomialMat.at(n - i, n - j);
            uint J = element.position({j}, n);
            Matrix.at(I, J) = binom * Bmoment(element.position({i + j}, n));
        }
    }
}