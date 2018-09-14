#include "MassM.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n) ((n + 1) * (n + 1)) // (((n+1) * (n+2)) / 2)
#ifdef LENB
#undef LENB
#endif
#define LENB(n) (2 * n + 2)

BMass2DTri::BMass2DTri(int q, int n)
    : BMoment2DTri(q, 2 * n),
      Matrix(LEN(n), LEN(n), arma::fill::none),
      BinomialMat(LENB(n), LENB(n), arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenMass = ((n + 2) * (n + 1)) / 2;
    lenBinomialMat = 2 * n + 2;
}

BMass2DTri::BMass2DTri(int q, int n, double T[][2])
    : BMoment2DTri(q, 2 * n, T),
      Matrix(LEN(n), LEN(n), arma::fill::none),
      BinomialMat(LENB(n), LENB(n), arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    lenMass = ((n + 2) * (n + 1)) / 2;
    lenBinomialMat = 2 * n + 2;

    setTriangle(T[0], T[1], T[2]);
}

BMass2DTri::~BMass2DTri()
{
}

void BMass2DTri::compute_binomials()
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

// this method was based on the code by M. Ainsworth
// G. Andriamro and O. Davydov
void BMass2DTri::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1.0 / BinomialMat.at(n, n);

    int M = MAX(2 * n, q - 1);

    int iMu = 0;
    for (int mu0 = n; mu0 >= 0; mu0--)
    {
        int iEta = 0;
        for (int eta0 = n; eta0 >= 0; eta0--)
        {
            int iMuEta0 = BMoment2DTri::position(mu0, eta0, M); // the Bmoments are indexed w.r.t position2d2(. , MAX(2n,q-1) )
            double w1 = Const * BinomialMat.at(mu0, eta0);

            int mu2 = 0;
            for (int mu1 = n - mu0; mu1 >= 0; mu1--, mu2++, iMu++)
            {
                int eta2 = 0;
                for (int eta1 = n - eta0; eta1 >= 0; eta1--, eta2++, iEta++)
                {
                    double w2 = w1 * BinomialMat.at(mu1, eta1) * BinomialMat.at(mu2, eta2);

                    int iMuEta = iMuEta0 + mu1 + eta1;

                    Matrix.at(iMu, iEta) = w2 * get_bmoment(iMuEta);
                }
                if (mu1 > 0)
                    iEta -= eta2;
            }
            if (eta0 > 0)
                iMu -= mu2;
        }
    }
}

/* void BMass2DTri::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1.0 / BinomialMat(n, n);

    int M = MAX(2 * n, q - 1); // constant for indexing

    for (int a1 = 0; a1 <= n; a1++) {
        for (int b1 = 0; b1 <= n; b1++) {
            double w1 = Const * BinomialMat(a1, b1);
            for (int a2 = 0; a2 <= n - a1; a2++){
                for (int b2 = 0; b2 <= n - b1; b2++) {
                    double w2 = w1 * BinomialMat(a2, b2);
                    int i = BMoment2DTri::position(a1, a2, n);
                    int j = BMoment2DTri::position(b1, b2, n);
                    int MomI = BMoment2DTri::position(a1 + b1, a2 + b2, M);
                    Matrix(i, j) = w2 * get_bmoment(MomI); // TODO: indexes missing
                } 
            }
        }
    }


} */

void BMass2DTri::compute_matrix(std::function<double(double, double)> f)
{
    setFunction(f);
    compute_matrix();
}

void BMass2DTri::compute_matrix(const arma::vec &Fval)
{
    setFunction(Fval);
    compute_matrix();
}