#include "Derivatives.h"
#include "JacobiGaussNodes.h"

#ifndef LEN
#define LEN(N) ((MAX(n, q - 1) + 1) * (MAX(n, q - 1) + 1))
#endif

using namespace TriD;
using namespace arma;

TriangleDerivative::TriangleDerivative(int q, int n) 
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

void TriangleDerivative::setTriangle(const mat &vertices)
{
    this->vertices = vertices;
}

void TriangleDerivative::setTriangle(double v1[], double v2[], double v3[])
{
    vertices(0, 0) = v1[0];
    vertices(0, 1) = v1[1];
    vertices(1, 0) = v2[0];
    vertices(1, 1) = v2[1];
    vertices(2, 0) = v3[0];
    vertices(2, 1) = v3[1];
}

void TriangleDerivative::setFunction(const mat &Fval)
{

}

void TriangleDerivative::compute_binomials(Mat<int64_t> &BinomialMat, int lenBinom)
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
    : TriangleDerivative(q, n) { }

void dXi_dXi::compute_matrix()
{
    int m = MAX(n+1, q);
    // calculate mu_0 = bmoment
    BMoment2DTri mu_0(q, 2 * (n - 1));
    mu_0.setFunction(Fval);
    mat triangle = {{0, 0}, {1, 0}, {0, 1}}; // TODO: change this triangle
    mu_0.setTriangle(triangle);
    // calculate mu_1
    vec mu_1(m * m);
    vec mu_1_inter(m * m);
    // calculate just like in BMoment2DTri, except this is a different form
    // convert first index
    // calculate mu_2
    // in O(n q^2 * n^2 q) each
    vec mu_2(m * m);
    vec mu_2_inter(m * m);

    double Const = n * n * (1.0 / BinomialMat(n-1, n-1));
    // compute matrix using mu_0, mu_1 and mu_2 in O(n^4)
}

/*****************************************
 ********** Defining dEta_dEta ***********
 *****************************************/

dEta_dEta::dEta_dEta(int q, int n)
    : TriangleDerivative(q, n) { }

void dEta_dEta::compute_matrix()
{
}

/****************************************
 ********** Defining dXi_dEta ***********
 ****************************************/

dXi_dEta::dXi_dEta(int q, int n)
    : TriangleDerivative(q, n) { }

void dXi_dEta::compute_matrix()
{

} 