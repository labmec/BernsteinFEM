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
        BinomialMat.at(i, 0) = 1;

    for (int j = 1; j < lenBinomialMat; j++)
        BinomialMat.at(0, j) = 1;

    for (int k = 1; k < lenBinomialMat; k++)
    {
        for (int l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat.at(k, l) += BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

void BStiff1D::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1. / BinomialMat.at(n - 1, n - 1) / pow(b - a, n);
    
    for (int i = 0; i < lenStiff - 1; i++)
    {
        for (int j = 0; j < lenStiff - 1; j++)
        {
            double w = BinomialMat.at(i, j) * Const;
            w *= BinomialMat.at(n - i - 1, n - j - 1);

            for (int k = 1; k <= 2; k++)
            {
                for (int l = 1; l <= 2; l++)
                {
                    int I = position(i + 2 - k, n);
                    int J = position(j + 2 - l, n);
                    Matrix.at(I, J) += (n * n) * w * grad(k, l) * get_bmoment(i + j);
                }
            }
        }
    }
}

void BStiff1D::compute_matrix(const arma::vec &Fval)
{
    setFunction(Fval);
    compute_matrix();
}

void BStiff1D::compute_matrix(std::function<double (double)> f)
{
    setFunction(f);
    compute_matrix();
}