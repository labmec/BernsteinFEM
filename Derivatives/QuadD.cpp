#include "Derivatives.h"

#ifndef LEN
#define LEN(N) ((n + 1) * (n + 1))
#endif

using namespace QuadD;
using namespace arma;

/****************************************
 ********** Defining dXi_dXi ************
 ****************************************/

void dXi_dXi::compute_matrix()
{
}

/****************************************
 ********** Defining dEta_dEta ************
 ****************************************/

void dEta_dEta::compute_matrix()
{
}

/****************************************
 ********** Defining dXi_dEta ************
 ****************************************/

void dXi_dEta::dXi_dEta(int q, int n) :
    Matrix(LEN(n), LEN(n), fill::none),
    BinomialMat(n + 1, n + 1, fill::zeros)
{
    this->q = q;
    this->n = n;

    len = LEN(n);
    lenBinom = n + 1;
    
    compute_binomials();
}

void dXi_dEta::compute_matrix()
{
    compute_moments();

    Matrix.zeros();

    // then arrange the terms
    double Const = n * n * (1.0 / BinomialMat(n, n));

    for (int a1 = 0; a1 < n; a1++)
    {
        for (int b1 = 0; b1 < n; b1++)
        {
            double w1 = BinomialMat(a1, b1);
            double _w1 = BinomialMat(a1 + 1, b1);
            for (int a2 = 0; a2 < n; a2++)
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

                    double w3 = BinomialMat
                }
            }
        }
    }
}