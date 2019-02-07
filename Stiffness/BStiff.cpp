#include "StiffM.h"

// default constructor
BStiff::BStiff(uint q, uint n)
    : Matrix(), BinomialMat(n + 1, n + 1, arma::fill::ones)
{
    this->q = q;
    this->n = n;
    lenBinomialMat = n + 1;
    computeBinomials();
}

// copy constructor
BStiff::BStiff(const BStiff &cp)
{
    q = cp.q;
    n = cp.n;
    lenBinomialMat = cp.lenBinomialMat;
    BinomialMat = cp.BinomialMat;
    lenStiff = cp.lenStiff;
    Matrix = cp.Matrix;
}

// copy assignment operator (unnecessary)
BStiff &BStiff::operator=(const BStiff &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        lenBinomialMat = cp.lenBinomialMat;
        BinomialMat = cp.BinomialMat;
        lenStiff = cp.lenStiff;
        Matrix = cp.Matrix;
    }
    return *this;
}

void BStiff::computeBinomials()
{
    for (uint k = 1; k < lenBinomialMat; k++)
    {
        for (uint l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat.at(k, l) = BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

// returns the length of the matrix
uint BStiff::length()
{
    return lenStiff;
}

// returns the value of Matrix[i][j]
double BStiff::getMatrixValue(uint i, uint j)
{
    return Matrix(i, j);
}

// returns the matrix
const arma::mat &BStiff::getMatrix()
{
    return Matrix;
}

// other methods

// assign 0 to all elements of the matrix
void BStiff::zero()
{
    Matrix.zeros();
}
