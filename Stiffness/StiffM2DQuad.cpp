#include "StiffM.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n) ((n + 1) * (n + 1))

BStiff2DQuad::BStiff2DQuad(int q, int n)
    : BMoment2DQuad(q, 2 * (n - 1)),
      Matrix(LEN(n), LEN(n), arma::fill::none),
      BinomialMat(n + 1, n + 1, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenStiff = (n + 1) * (n + 1);
    lenBinomialMat = n + 1;
}

BStiff2DQuad::~BStiff2DQuad()
{
}

void BStiff2DQuad::compute_binomials()
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

void BStiff2DQuad::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1 / (BinomialMat(n - 1, n - 1) * BinomialMat(n - 1, n - 1));

    for (int i = 0; i < lenStiff - 1; i++)
    {
        for (int j = 0; j < lenStiff - 1; j++)
        {
            double w = (BinomialMat(i, j) * BinomialMat(i, j)) * Const;
            w *= (BinomialMat(n - i - 1, n - j - 1) * BinomialMat(n - i - 1, n - j - 1));

            for (int k = 1; k <= 2; k++)
                for (int l = 1; l <= 2; l++)
                    for (int k2 = 1; k2 <= 2; k2++)
                        for (int l2 = 1; l2 <= 2; l2++)
                        {
                            double grad = BStiff1D::grad(k, l) * BStiff1D::grad(k2, l2);
                            Matrix(i, j) += (n * n) * w * grad * get_bmoment(i + j);
                        }
        }
    }
}

void BStiff2DQuad::compute_matrix(arma::vec Fval)
{
    setFunction(Fval);
    compute_matrix();
}

void BStiff2DQuad::compute_matrix(double (*f)(double, double))
{
    setFunction(f);
    compute_matrix();
}