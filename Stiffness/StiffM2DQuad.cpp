#include "StiffM.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n) ((n + 1) * (n + 1))

#ifndef USE_DERIVATIVES

BStiff2DQuad::BStiff2DQuad(int q, int n)
    : Matrix(LEN(n), LEN(n), arma::fill::none),
      Moments1(q, 2 * (n - 1), 2 * n),
      Moments2(q, 2 * n, 2 * (n - 1)),
      BinomialMat(n + 1, n + 1, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenStiff = LEN(n);
    lenBinomialMat = n + 1;
}

BStiff2DQuad::~BStiff2DQuad()
{
}

void BStiff2DQuad::setFunction(arma::mat Fval) {
    Moments1.setFunction(Fval);
    Moments2.setFunction(Fval);
}

void BStiff2DQuad::setFunction(arma::vec Fval) {
    Moments1.setFunction(Fval);
    Moments2.setFunction(Fval);
}

void BStiff2DQuad::setFunction(std::function<double (double, double)> f) {
    Moments1.setFunction(f) ;
    Moments2.setFunction(f) ;
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

// TODO: update this function;
// either make it work this way or change it to compute using the QuadD namespace
void BStiff2DQuad::compute_matrix()
{
    Moments1.compute_moments();
    Moments2.compute_moments();
    arma::mat m1 = Moments1.get_bmoment();
    arma::mat m2 = Moments2.get_bmoment();
    compute_binomials();
    
    double Const = n * n * (1.0  / (BinomialMat(n - 1, n - 1) * BinomialMat(n - 1, n - 1)));

    for (int a1 = 0; a1 < n; a1++)
    {
        for (int b1 = 0; b1 < n; b1++)
        {
            double w1 = Const * BinomialMat(a1, b1) * BinomialMat(n - a1, n - b1);
            for (int a2 = 0; a2 <= n; a2++)
            {
                for (int b2 = 0; b2 <= n; b2++)
                {
                    double w2 = w1 * BinomialMat(a2, b2) * BinomialMat(n - a2, n - b2);
                    double mom1 = w2 * m1( BMoment2DQuad::position(a1 + b1, a2 + b2, n) );
                    double mom2 = w2 * m2( BMoment2DQuad::position(a1 + b1, a2 + b2, n) );
                    double mom_tot = mom1 + mom2;

                    int i = BMoment2DQuad::position(a1, a2, n);
                    int j = BMoment2DQuad::position(b1, b2, n);
                    int I = BMoment2DQuad::position(a1 + 1, a2, n);
                    int J = BMoment2DQuad::position(b1 + 1, b2, n);

                    Matrix(i, j) += mom_tot;

                    Matrix(I, J) += mom_tot;

                    Matrix(I, j) += mom_tot;

                    Matrix(i, J) += mom_tot;
                }
            }
        }
    }
}

void BStiff2DQuad::compute_matrix(arma::vec Fval)
{
    setFunction(Fval);
    compute_matrix();
}

void BStiff2DQuad::compute_matrix(std::function<double (double, double)> f)
{
    setFunction(f);
    compute_matrix();
}

#endif