#include "StiffM.h"

BStiff1D::BStiff1D(int q, int n) : BMoment1D(q, 2 * (n - 1))
{
    this->q = q;
    this->n = n;

    lenStiff = n + 1;
    lenBinomialMat = n + 1;

    Matrix = create_matrix();
    BinomialMat = create_binomialMat();
}

BStiff1D::~BStiff1D()
{
    delete_matrix(Matrix);
    delete_binomialMat(BinomialMat);
}

double **BStiff1D::create_matrix()
{
    double *aux = new double[lenStiff * lenStiff];
    double **matrix = new double *[lenStiff];

    for (int i = 0; i < lenStiff; aux += lenStiff, i++)
        matrix[i] = aux;

    return matrix;
}

void BStiff1D::delete_matrix(double **matrix)
{
    delete matrix[0];
    delete matrix;
}

int **BStiff1D::create_binomialMat()
{
    int *aux = new int[lenBinomialMat * lenBinomialMat];
    int **matrix = new int *[lenBinomialMat];

    for (int i = 0; i < lenBinomialMat; aux += lenBinomialMat, i++)
        matrix[i] = aux;

    return matrix;
}

void BStiff1D::delete_binomialMat(int **binomialMat)
{
    delete binomialMat[0];
    delete binomialMat;
}

void BStiff1D::compute_binomials()
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

void BStiff1D::compute_matrix()
{
    compute_moments();
    compute_binomials();

    double Const = 1. / BinomialMat[n - 1][n - 1];
    
    for (int i = 0; i < lenStiff - 1; i++)
    {
        for (int j = 0; j < lenStiff - 1; j++)
        {
            double w = BinomialMat[i][j] * Const;
            w *= BinomialMat[n - i - 1][n - j - 1];

            for (int k = 1; k <= 2; k++)
                for (int l = 1; l <= 2; l++)
                    Matrix[i + 2 - k][j + 2 - l] += (n * n) * w * grad(k, l) * get_bmoment(i + j);
        }
    }
}

void BStiff1D::compute_matrix(double *Fval)
{
    setFunction(Fval);
    compute_matrix();
}

void BStiff1D::compute_matrix(double (*f)(double))
{
    setFunction(f);
    compute_matrix();
}