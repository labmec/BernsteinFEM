#include "MassM.h"

BMass::BMass(uint q, uint n)
    : Matrix(), BinomialMat(n + 1, n + 1, arma::fill::ones)
{
    this->q = q;
    this->n = n;
    lenBinomialMat = n + 1;
    computeBinomials();
}

BMass::BMass(const BMass &cp)
{
    q = cp.q;
    n = cp.n;
    lenMass = cp.lenMass;
    lenBinomialMat = cp.lenBinomialMat;
    Matrix = cp.Matrix;
    BinomialMat = cp.BinomialMat;
}

BMass &BMass::operator=(const BMass &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        lenMass = cp.lenMass;
        lenBinomialMat = cp.lenBinomialMat;
        Matrix = cp.Matrix;
        BinomialMat = cp.BinomialMat;
    }
    return *this;
}

void BMass::computeBinomials()
{
    for (uint k = 1; k < lenBinomialMat; k++)
    {
        for (uint l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat.at(k, l) = BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

// getters

uint BMass::length()
{
    return lenMass;
}

// returns the polynomial order of the basis
uint BMass::getPOrder()
{
    return n;
}

// returns a reference to the Mass Matrix
const arma::mat &BMass::getMatrix()
{
    return Matrix;
}

// returns the value of Matrix[i][j]
double BMass::getMatrixValue(uint i, uint j)
{
    try
    {
        return Matrix(i, j);
    }
    catch (std::logic_error &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Max value for arguments: " << Matrix.n_cols << std::endl;
        throw std::logic_error(e);
    }
};