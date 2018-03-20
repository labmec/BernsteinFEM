#include "MassM.h"

BMass2DTri::BMass2DTri(int q, int n) : BMoment2DTri(q, 2 * n)
{
    this->q = q;
    this->n = n;

    lenMass = ((n + 2) * (n + 1)) / 2;
    lenBinomialMat = 2 * n + 2;

    Matrix = create_matrix();
    BinomialMat = create_binomialMat();
}

BMass2DTri::BMass2DTri(int q, int n, double T[][2]) : BMoment2DTri(q, 2 * n, T)
{
    this->q = q;
    this->n = n;

    lenMass = ((n + 2) * (n + 1)) / 2;
    lenBinomialMat = 2 * n + 2;

    Matrix = create_matrix();
    BinomialMat = create_binomialMat();
}

BMass2DTri::~BMass2DTri()
{
    delete_matrix(Matrix);
    delete_binomialMat(BinomialMat);
}

double **BMass2DTri::create_matrix()
{
    double *aux = new double[lenMass * lenMass];
    double **matrix = new double *[lenMass];

    for (int i = 0; i < lenMass; aux += lenMass, i++)
        matrix[i] = aux;

    return matrix;
}

void BMass2DTri::delete_matrix(double **matrix)
{
    delete matrix[0];
    delete matrix;
}

int **BMass2DTri::create_binomialMat()
{
    int *aux = new int[lenBinomialMat * lenBinomialMat];
    int **matrix = new int *[lenBinomialMat];

    for (int i = 0; i < lenBinomialMat; aux += lenBinomialMat, i++)
        matrix[i] = aux;

    return matrix;
}

void BMass2DTri::delete_binomialMat(int **binomialMat)
{
    delete binomialMat[0];
    delete binomialMat;
}

void BMass2DTri::compute_binomials()
{
    for (int i = 0; i < lenBinomialMat; i++)
    {
        for (int j = 0; j < lenBinomialMat; j++)
        {
            BinomialMat[i][j] = 0;
        }
    }
    for (int i = 0; i < lenBinomialMat; i++)
        BinomialMat[i][0] += 1;

    for (int j = 1; j < lenBinomialMat; j++)
        BinomialMat[0][j] += 1;

    for (int k = 1; k < lenBinomialMat; k++)
    {
        for (int l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat[k][l] += BinomialMat[k][l - 1] + BinomialMat[k - 1][l];
        }
    }
}

void BMass2DTri::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1.0 / BinomialMat[n][n];

    // initialize the matrix with 0's
    for (int i = 0; i < lenMass; i++)
    {
        for (int j = 0; j < lenMass; j++)
            Matrix[i][j] = 0;
    }

    int M = MAX(2 * n, q - 1);

    int I = 0;
    for (int mu0 = n; mu0 >= 0; mu0--)
    {
        int J = 0;
        for (int eta0 = n; eta0 >= 0; eta0--)
        {
            int iMuEta0 = (mu0 + eta0) * (M + 1); // the Bmoments are indexed w.r.t position2d2(. , MAX(2n,q-1) )
            double w1 = Const * BinomialMat[mu0][eta0];

            int mu2 = 0;
            for (int mu1 = n - mu0; mu1 >= 0; mu1--, mu2++, I++)
            {
                int eta2 = 0;
                for (int eta1 = n - eta0; eta1 >= 0; eta1--, eta2++, J++)
                {
                    double w2 = w1 * BinomialMat[mu1][eta1] * BinomialMat[mu2][eta2];

                    int iMuEta = iMuEta0 + mu1 + eta1;

                    Matrix[I][J] = w2 * get_bmoment(iMuEta);
                }
                if (mu1 > 0)
                    J -= eta2;
            }
            if (eta0 > 0)
                I -= mu2;
        }
    }
}

void BMass2DTri::compute_matrix(double (*f)(double, double))
{
    setFunctionDef(f);
    compute_matrix();
}

void BMass2DTri::compute_matrix(double *Fval)
{
    setFunctionValue(Fval);
    compute_matrix();
}