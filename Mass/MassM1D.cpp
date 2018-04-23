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

void BMass1D::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1.0 / BinomialMat(n, n);
    double binom;

    for (int i = 0; i < lenMass; i++)
    {
        for (int j = 0; j < lenMass; j++)
        {
            binom = Const * BinomialMat(i, j) * BinomialMat(n - i, n - j);
            Matrix(i, j) = binom * get_bmoment(i + j);
        }
    }
}

void BMass1D::compute_matrix(double (*f)(double))
{
    setFunction(f);
    compute_matrix();
}

void BMass1D::compute_matrix(arma::vec Fval)
{
    setFunction(Fval);
    compute_matrix();
}