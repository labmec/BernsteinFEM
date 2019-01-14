#include "StiffM.h"

// default constructor
BStiff::BStiff(int q, int n)
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

inline
void BStiff::computeBinomials()
{
    for (int k = 1; k < lenBinomialMat; k++)
    {
        for (int l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat.at(k, l) = BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

// returns the length of the matrix
inline
int BStiff::length()
{
    return lenStiff;
}

// returns the value of Matrix[i][j]
inline
double BStiff::getMatrixValue(int i, int j)
{
    return Matrix(i, j);
}

// returns the matrix
inline
const arma::mat &BStiff::getMatrix()
{
    return Matrix;
}

// other methods

// assign 0 to all elements of the matrix
inline
void BStiff::zero()
{
    Matrix.zeros();
}
