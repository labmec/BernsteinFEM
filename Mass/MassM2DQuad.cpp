#include "MassM.h"

BMass2DQuad::BMass2DQuad(int q, int n)
    : BMoment2DQuad(q, 2 * n),
      Matrix((n + 1) * (n + 1), (n + 1) * (n + 1), arma::fill::none),
      BinomialMat(n * 2 + 2, n * 2 + 2, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenMass = (n + 1) * (n + 1);
    lenBinomialMat = 2 * n + 2;
}

BMass2DQuad::~BMass2DQuad()
{
}

void BMass2DQuad::compute_binomials()
{
    for (int i = 0; i < lenBinomialMat; i++)
        BinomialMat(i, 0) += 1;

    for (int j = 1; j < lenBinomialMat; j++)
        BinomialMat(0, j) += 1;

    for (int k = 1; k < lenBinomialMat; k++)
    {
        for (int l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat(k, l) += BinomialMat(k, l - 1) + BinomialMat(k - 1, l);
        }
    }
}

void BMass2DQuad::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1.0 / (BinomialMat(n, n) * BinomialMat(n, n));

    // since it is a simple tensor product, this is just like in the 1D case
    // except it's doubled
    for (int i = 0; i < lenMass; i++)
    {
        for (int j = 0; j < lenMass; j++)
        {
            double w = Const * BinomialMat(i, j) * BinomialMat(i, j);
            w *= (BinomialMat(n - i, n - j) * BinomialMat(n - i, n - j));
            Matrix(i, j) = w * get_bmoment(i + j);
        }
    }
}

void BMass2DQuad::compute_matrix(double (*f)(double, double))
{
    setFunction(f);
    compute_matrix();
}

void BMass2DQuad::compute_matrix(arma::vec Fval)
{
    setFunction(Fval);
    compute_matrix();
}