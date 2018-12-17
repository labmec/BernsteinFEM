#include <iostream>
using std::cin;
using std::cout;
using std::endl;

#include "MassM.h"

// maybe trade q to 2*q in the base class constructor
BMass1D::BMass1D(int q, int n, Element<Element_t::LinearEl> el)
    : BMass(q, n), BMoment1D(q, 2 * n, el)
{
    lenMass = n + 1;
    Matrix.set_size(lenMass, lenMass);
}

BMass1D::~BMass1D()
{
}

void BMass1D::computeMatrix()
{
    int n = BMass::n;
    computeMoments();

    double Const = 1.0 / BinomialMat.at(n, n);

    for (int i = 0; i < lenMass; i++)
    {
        int I = position(i, n);
        for (int j = 0; j < lenMass; j++)
        {
            double binom = Const * BinomialMat.at(i, j) * BinomialMat.at(n - i, n - j);
            int J = position(j, n);
            Matrix.at(I, J) = binom * Bmoment(position(i + j, n));
        }
    }
}