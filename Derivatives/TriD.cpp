#include "Derivatives.h"
#include "JacobiGaussNodes.h"

#ifndef LEN
#define LEN(N) ((MAX(n, q - 1) + 1) * (MAX(n, q - 1) + 1))
#endif

using namespace TriD;
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
    lenBinom = n + 1; // TODO: verificar esse valor

    compute_binomials(BinomialMat, lenBinom);
}

void dXi_dXi::compute_matrix()
{
}

/*****************************************
 ********** Defining dEta_dEta ***********
 *****************************************/

dEta_dEta::dEta_dEta(int q, int n)
    : Matrix(LEN(n), LEN(n), fill::none),
      Fval(q, q, fill::none),
      BinomialMat(n + 1, n + 1, fill::zeros)
{
    this->q = q;
    this->n = n;

    len = LEN(n);
    lenBinom = n + 1; // TODO: verificar esse valor

    compute_binomials(BinomialMat, lenBinom);
}

void dEta_dEta::compute_matrix()
{
}

/****************************************
 ********** Defining dXi_dEta ***********
 ****************************************/

dXi_dEta::dXi_dEta(int q, int n)
    : Matrix(LEN(n), LEN(n), fill::none),
      Fval(q, q, fill::none),
      BinomialMat(n + 1, n + 1, fill::zeros)
{
    this->q = q;
    this->n = n;

    len = LEN(n);
    lenBinom = n + 1; // TODO: verificar esse valor

    compute_binomials(BinomialMat, lenBinom);
}

void dXi_dEta::compute_matrix()
{

} 