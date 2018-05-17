#include "Derivatives.h"
#include "JacobiGaussNodes.h"

#include <iostream>

#ifndef LEN
#define LEN(N) ((N + 1) * (N + 1))
#endif

using namespace QuadD;
using namespace arma;

void compute_binomials(Mat<int64_t> &BinomialMat, int lenBinom)
{
    for (int i = 0; i < lenBinom; i++)
        BinomialMat(i, 0) = 1;

    for (int j = 1; j < lenBinom; j++)
        BinomialMat(0, j) = 1;

    for (int k = 1; k < lenBinom; k++)
    {
        for (int l = 1; l < lenBinom; l++)
        {
            BinomialMat(k, l) += BinomialMat(k, l - 1) + BinomialMat(k - 1, l);
        }
    }
}

/****************************************
 ********** Defining dXi_dXi ************
 ****************************************/

dXi_dXi::dXi_dXi(int q, int n)
    : Matrix(LEN(n), LEN(n), fill::none),
      Fval(q, q, fill::none),
      BinomialMat(n + 1, n + 1, fill::zeros)
{
    this->q = q;
    this->n = n;

    len = LEN(n);
    lenBinom = n + 1;

    compute_binomials(BinomialMat, lenBinom);
}

void dXi_dXi::setFunction(Mat<double> FVal)
{
    for (int i = 0; i < q; i++)
        for (int j = 0; j < q; j++)
            Fval(i, j) = FVal(i, j);
}

void dXi_dXi::compute_matrix()
{
    Matrix.zeros();

    // first we compute the integrals
    int nXi = 2 * (n - 1);
    int nEta = 2 * n;
    Mat<double> moments(nXi + 1, nEta + 1, fill::zeros);

    for (int iXi = 0; iXi < q; iXi++)
    {
        double xi = (legendre_xi(q, iXi) + 1.0) * 0.5; // the master element is [0, 1] X [0, 1]
        double sXi = 1.0 - xi;
        double rXi = xi / sXi;
        double omegaXi = legendre_w(q, iXi) * 0.5;

        double wXi = omegaXi * pow(sXi, nXi); // O(q * n)
        for (int aXi = 0; aXi <= nXi; aXi++)
        {
            for (int iEta = 0; iEta < q; iEta++)
            {
                double eta = (legendre_xi(q, iEta) + 1.0) * 0.5;
                double sEta = 1.0 - eta;
                double rEta = eta / sEta;
                double omegaEta = legendre_w(q, iEta) * 0.5;

                double wEta = omegaEta * pow(sEta, nEta); // O(q^2 * n^2)
                for (int aEta = 0; aEta <= nEta; aEta++)
                {
                    double w = wXi * wEta * Fval(iXi, iEta);
                    moments(aXi, aEta) += w;

                    wEta *= rEta * (nEta - aEta) / (1.0 + aEta);
                }
            }
            wXi *= rXi * (nXi - aXi) / (1.0 + aXi);
        }
    } // O(q^2 * n^2)

    // test if the moments were calculated correctly
    /* for (int i = 0; i <= nXi; i++) {
        for(int j = 0; j <= nEta; j++)
        {
            std::cout << std::scientific << moments(i, j) << ",\t";
        }
        std::cout << std::endl;
    } */

    // then we arrange the terms

    double Const = n * n * (1.0 / (BinomialMat(n, n) * BinomialMat(n - 1, n - 1)));

    for (int a1 = 0; a1 < n; a1++)
    {
        for (int b1 = 0; b1 < n; b1++)
        {
            double w1 = Const * BinomialMat(a1, b1);
            for (int a2 = 0; a2 <= n; a2++)
            {
                for (int b2 = 0; b2 <= n; b2++)
                {
                    double w2 = BinomialMat(a2, b2);
                    double mom = moments(a1 + b1, a2 + b2);

                    int i = BMoment2DQuad::position(a1, a2, n);
                    int j = BMoment2DQuad::position(b1, b2, n);
                    int I = BMoment2DQuad::position(a1 + 1, a2, n);
                    int J = BMoment2DQuad::position(b1 + 1, b2, n);

                    Matrix(i, j) += mom * w1 * w2;

                    Matrix(I, J) += mom * w1 * w2;

                    Matrix(I, j) -= mom * w1 * w2;

                    Matrix(i, J) -= mom * w1 * w2;
                }
            }
        }
    } // O(n^4)
}

/****************************************
 ********** Defining dEta_dEta ************
 ****************************************/

dEta_dEta::dEta_dEta(int q, int n)
    : Matrix(LEN(n), LEN(n), fill::none),
      Fval(q, q, fill::none),
      BinomialMat(n + 1, n + 1, fill::zeros)
{
    this->q = q;
    this->n = n;

    len = LEN(n);
    lenBinom = n + 1;

    compute_binomials(BinomialMat, lenBinom);
}

void dEta_dEta::setFunction(Mat<double> FVal)
{
    for (int i = 0; i < q; i++)
        for (int j = 0; j < q; j++)
            Fval(i, j) = FVal(i, j);
}

void dEta_dEta::compute_matrix()
{
    // dEta_dEta is basically the same as dXi_dXi, we could do it by transposing the FVal matrix
    // and then transposing the resulting matrix
    Matrix.zeros();

    // first we compute the integrals
    int nXi = 2 * n;
    int nEta = 2 * (n - 1);
    Mat<double> moments(nXi + 1, nEta + 1, fill::zeros);

    for (int iXi = 0; iXi < q; iXi++)
    {
        double xi = (legendre_xi(q, iXi) + 1.0) * 0.5; // the master element is [0, 1] X [0, 1]
        double sXi = 1 - xi;
        double rXi = xi / sXi;
        double omegaXi = legendre_w(q, iXi) * 0.5;

        double wXi = omegaXi * pow(sXi, nXi);
        for (int aXi = 0; aXi <= nXi; aXi++)
        {
            for (int iEta = 0; iEta < q; iEta++)
            {
                double eta = (legendre_xi(q, iEta) + 1.0) * 0.5;
                double sEta = 1 - eta;
                double rEta = eta / sEta;
                double omegaEta = legendre_w(q, iEta) * 0.5;

                double wEta = omegaEta * pow(sEta, nEta);
                for (int aEta = 0; aEta <= nEta; aEta++)
                {
                    moments(aXi, aEta) += wXi * wEta * Fval(iXi, iEta);

                    wEta *= rEta * ((nEta - aEta) / (1 + aEta));
                }
            }
            wXi *= rXi * ((nXi - aXi) / (1 + aXi));
        }
    } // O(q^2 * n^2)

    // then we arrange the terms

    double Const = n * n * (1.0 / (BinomialMat(n, n) * BinomialMat(n - 1, n - 1)));

    for (int a1 = 0; a1 <= n; a1++)
    {
        for (int b1 = 0; b1 < n; b1++)
        {
            double w1 = BinomialMat(a1, b1);
            for (int a2 = 0; a2 <= n; a2++)
            {
                for (int b2 = 0; b2 < n; b2++)
                {
                    double w2 = Const * BinomialMat(a2, b2);
                    double mom = moments(a1 + b1, a2 + b2);

                    int i = BMoment2DQuad::position(a1, a2, n);
                    int j = BMoment2DQuad::position(b1, b2, n);
                    int I = BMoment2DQuad::position(a1, a2 + 1, n);
                    int J = BMoment2DQuad::position(b1, b2 + 1, n);

                    Matrix(i, j) += mom * w1 * w2;

                    w1 = BinomialMat(a1 + 1, b1 + 1);
                    Matrix(I, J) += mom * w1 * w2;

                    w1 = BinomialMat(a1 + 1, b1);
                    Matrix(I, j) -= mom * w1 * w2;

                    w1 = BinomialMat(a1, b1 + 1);
                    Matrix(i, J) -= mom * w1 * w2;
                }
            }
        }
    } // O(n^4)
}

/****************************************
 ********** Defining dXi_dEta ************
 ****************************************/

dXi_dEta::dXi_dEta(int q, int n)
    : BMoment2DQuad(q, 2 * n - 1),
      Matrix(LEN(n), LEN(n), fill::none),
      BinomialMat(n + 1, n + 1, fill::zeros)
{
    this->q = q;
    this->n = n;

    len = LEN(n);
    lenBinom = n + 1;

    compute_binomials(BinomialMat, lenBinom);
}

void dXi_dEta::compute_matrix()
{
    compute_moments();
    Matrix.zeros();

    // then arrange the terms
    double Const = n * n * (1.0 / BinomialMat(n, n));

    for (int a1 = 0; a1 < n; a1++)
    {
        for (int b1 = 0; b1 <= n; b1++)
        {
            double w1 = BinomialMat(a1, b1);
            for (int a2 = 0; a2 <= n; a2++)
            {
                for (int b2 = 0; b2 < n; b2++)
                {
                    double w2 = BinomialMat(a2, b2);
                    double mom = get_bmoment(a1 + b1, a2 + b2);

                    int i = position(a1, a2, n);
                    int j = position(b1, b2, n);
                    int I = position(a1 + 1, a2, n); // a1 + 1
                    int J = position(b1, b2 + 1, n); // b2 + 1

                    Matrix(i, j) += Const * w1 * w2 * mom;
                    Matrix(I, J) += Const * w1 * w2 * mom;

                    Matrix(I, j) -= Const * w1 * w2 * mom;
                    Matrix(i, J) -= Const * w1 * w2 * mom;

                    std::cout << "Matrix max: " << Matrix.n_cols << std::endl;
                    std::cout << "i:\t" << i << endl;
                    std::cout << "j:\t" << j << endl;
                    std::cout << "I:\t" << I << endl;
                    std::cout << "J:\t" << J << endl;
                }
            }
        }
    }
} // takes the time to compute moments with parameters q and 2 * n - 2 (q^2 * (2 * n - 2)^2) plus n^4
