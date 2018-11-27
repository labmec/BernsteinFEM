#include <iostream>
using std::cin;
using std::cout;
using std::endl;

#include "MassM.h"

// maybe trade q to 2*q in the base class constructor
BMass1D::BMass1D(int q, int n)
    : BMoment1D(q, 2 * n),
      Matrix(n + 1, n + 1, arma::fill::zeros),
      BinomialMat(n + 1, n + 1, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenMass = n + 1;
    lenBinomialMat = n + 1;
}

BMass1D::~BMass1D()
{
}

void BMass1D::compute_binomials()
{
    for (int i = 0; i < lenBinomialMat; i++)
        BinomialMat.at(i, 0) += 1;

    for (int j = 1; j < lenBinomialMat; j++)
        BinomialMat.at(0, j) += 1;

    for (int k = 1; k < lenBinomialMat; k++)
    {
        for (int l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat.at(k, l) += BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

void BMass1D::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1.0 / BinomialMat.at(n, n);

    for (int i = 0; i < lenMass; i++)
    {
        int I = position(i, n);
        for (int j = 0; j < lenMass; j++)
        {
            double binom = Const * BinomialMat.at(i, j) * BinomialMat.at(n - i, n - j);
            int J = position(j, n);
            Matrix.at(I, J) = binom * get_bmoment(position(i + j, n));
        }
    }
}

void BMass1D::compute_matrix(std::function<double(double)> f)
{
    setFunction(f);
    compute_matrix();
}

void BMass1D::compute_matrix(const arma::vec &Fval)
{
    setFunction(Fval);
    compute_matrix();
}