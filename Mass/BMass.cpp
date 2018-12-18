#include "MassM.h"

BMass::BMass(int q, int n)
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

inline
void BMass::computeBinomials()
{
    for (int k = 1; k < lenBinomialMat; k++)
    {
        for (int l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat.at(k, l) = BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

// getters

inline
int BMass::length()
{
    return lenMass;
}

// returns the polynomial order of the basis
inline
int BMass::getPOrder()
{
    return n;
}

// returns a reference to the Mass Matrix
inline
const arma::mat &BMass::getMatrix()
{
    return Matrix;
}

// returns the value of Matrix[i][j]
inline
double BMass::getMatrixValue(int i, int j)
{
    try
    {
        return Matrix(i, j);
    }
    catch (std::logic_error &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "Max value for arguments: " << Matrix.n_cols << std::endl;
    }
};