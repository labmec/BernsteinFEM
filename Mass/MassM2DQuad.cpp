#include "MassM.h"

#define LEN(n) ((n + 1) * (n + 1))

BMass2DQuad::BMass2DQuad(int q, int n)
    : BMoment2DQuad(q, 2 * n),
      Matrix(LEN(n), LEN(n), arma::fill::none),
      BinomialMat(n + 1, n + 1, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenMass = LEN(n);
    lenBinomialMat = n + 1;
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

    double Const = 1.0 / (BinomialMat(n, n) * BinomialMat(n, n)); // constant due to integration

    // since it is a simple tensor product, this is just like in the 1D case
    // except for indexing
    for (int a1 = 0; a1 <= n; a1++)
    {
        for (int a2 = 0; a2 <= n; a2++) 
        {
            for (int b1 = 0; b1 <= n; b1++)
            {
                for (int b2 = 0; b2 <= n; b2++)
                {
                    double w = Const * BinomialMat(a1, b1) * BinomialMat(a2, b2);
                    w *= (BinomialMat(n - a1, n - b1) * BinomialMat(n - a2, n - b2));
                    int i = position(a1, a2, n);
                    int j = position(b1, b2, n);
                    Matrix(i, j) = w * get_bmoment(a1 + b1, a2 + b2);
                }
            }
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