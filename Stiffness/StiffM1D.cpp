#include "StiffM.h"

BStiff1D::BStiff1D(int q, int n)
    : BMoment1D(q, 2 * (n - 1)),
      Matrix(n + 1, n + 1, arma::fill::none),
      BinomialMat(n + 1, n + 1, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenStiff = n + 1;
    lenBinomialMat = n + 1;
}

BStiff1D::~BStiff1D()
{
}

void BStiff1D::compute_binomials()
{
    for (int i = 0; i < lenBinomialMat; i++)
        BinomialMat(i, 0) = 1;

    for (int j = 1; j < lenBinomialMat; j++)
        BinomialMat(0, j) = 1;

    for (int k = 1; k < lenBinomialMat; k++)
    {
        for (int l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat(k, l) += BinomialMat(k, l - 1) + BinomialMat(k - 1, l);
        }
    }
}

void BStiff1D::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1. / BinomialMat(n - 1, n - 1) / pow(b - a, n);
    
    for (int i = 0; i < lenStiff - 1; i++)
    {
        for (int j = 0; j < lenStiff - 1; j++)
        {
            double w = BinomialMat(i, j) * Const;
            w *= BinomialMat(n - i - 1, n - j - 1);

            for (int k = 1; k <= 2; k++)
                for (int l = 1; l <= 2; l++)
                    Matrix(i + 2 - k, j + 2 - l) += (n * n) * w * grad(k, l) * get_bmoment(i + j);
        }
    }
}

void BStiff1D::compute_matrix(arma::vec Fval)
{
    setFunction(Fval);
    compute_matrix();
}

void BStiff1D::compute_matrix(std::function<double (double)> f)
{
    setFunction(f);
    compute_matrix();
}