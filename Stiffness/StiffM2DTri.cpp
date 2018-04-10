#include "StiffM.h"

BStiff2DTri::BStiff2DTri(int q, int n) : BMoment2DTri(q, 2 * (n - 1), 6)
{
    this->q = q;
    this->n = n;

    Moments = new BMoment2DTri(q, 2 * (n - 1), 3);

    lenStiff = (n + 1) * (n + 1);
    lenBinomialMat = n + 1;

    Matrix = create_matrix();
    BinomialMat = create_binomialMat();
}

BStiff2DTri::~BStiff2DTri()
{
    delete_matrix(Matrix);
    delete_binomialMat(BinomialMat);
    delete Moments;
}

double **BStiff2DTri::create_matrix()
{
    double *aux = new double[lenStiff * lenStiff];
    double **matrix = new double *[lenStiff];

    for (int i = 0; i < lenStiff; aux += lenStiff, i++)
        matrix[i] = aux;

    return matrix;
}

void BStiff2DTri::delete_matrix(double **matrix)
{
    delete matrix[0];
    delete matrix;
}

int **BStiff2DTri::create_binomialMat()
{
    int *aux = new int[lenBinomialMat * lenBinomialMat];
    int **matrix = new int *[lenBinomialMat];

    for (int i = 0; i < lenBinomialMat; aux += lenBinomialMat, i++)
        matrix[i] = aux;

    return matrix;
}

void BStiff2DTri::delete_binomialMat(int **binomialMat)
{
    delete binomialMat[0];
    delete binomialMat;
}

void BStiff2DTri::compute_binomials()
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

void BStiff2DTri::compute_normals()
{
    normalMat[0][0] = v2[1] - v3[1];
    normalMat[0][1] = v3[0] - v2[0];
    normalMat[1][0] = v3[1] - v1[1];
    normalMat[1][1] = v1[0] - v3[0];
    normalMat[2][0] = v1[1] - v2[1];
    normalMat[2][1] = v2[0] - v1[0];
}

void BStiff2DTri::compute_matrix()
{
    Moments->compute_moments();
    compute_normals();
    transform_BmomentC_Stiff2d(Moments, normalMat);

    // initialize matrix with 0's
    for (int i = 0; i < lenStiff; i++)
        for (int j = 0; j < lenStiff; j++)
            Matrix[i][j] = 0.0;

    double Const = n * n / 4. / Area2d(v1, v2, v3) / Area2d(v1, v2, v3) / BinomialMat[n - 1][n - 1]; // taking account of scaling between normals and gradients

    int m = n - 1;
    int M = MAX(2 * n - 2, q - 1);

    // the same loop as in Mass, but with m=n-1 instead of n, and with moving stencil
    // of position2ds in the stiffness matrix
    int iMu = 0;
    int iMu1 = 1;
    int iMu2 = 2;

    for (int mu0 = m; mu0 >= 0; mu0--, iMu1++, iMu2++)
    {
        int iEta = 0;
        int iEta1 = 1;
        int iEta2 = 2;

        for (int eta0 = m; eta0 >= 0; eta0--, iEta1++, iEta2++)
        {
            int iMuEta0 = (mu0 + eta0) * (M + 1);

            double v = BinomialMat[mu0][eta0]; //to reduce the number of accesses to binomialMat, for 3D and higher
            //this way the number of multiplications may be reduced, too

            int mu2 = 0; // mu1=m-mu0-mu2
            for (int mu1 = m - mu0; mu1 >= 0; mu1--, mu2++, iMu++, iMu1++, iMu2++)
            {

                int eta2 = 0; // eta2=m-eta0
                for (int eta1 = m - eta0; eta1 >= 0; eta1--, eta2++, iEta++, iEta1++, iEta2++)
                {

                    double w = v * BinomialMat[mu1][eta1] * BinomialMat[mu2][eta2]; //use Pascal directly

                    w *= Const;
                    int iMuEta;

                    iMuEta = iMuEta0 + mu1 + eta1;

                    Matrix[iMu][iEta] += w * get_bmoment(iMuEta, 0);  //alfa=[1,0,0], beta=[1,0,0]
                    Matrix[iMu][iEta1] += w * get_bmoment(iMuEta, 1); //beta=[0,1,0]
                    Matrix[iMu][iEta2] += w * get_bmoment(iMuEta, 2); //beta=[0,0,1]

                    Matrix[iMu1][iEta] += w * get_bmoment(iMuEta, 1);  //alfa=[0,1,0], beta=[1,0,0]
                    Matrix[iMu1][iEta1] += w * get_bmoment(iMuEta, 3); //beta=[0,1,0]
                    Matrix[iMu1][iEta2] += w * get_bmoment(iMuEta, 4); //beta=[0,0,1]

                    Matrix[iMu2][iEta] += w * get_bmoment(iMuEta, 2);  //alfa=[0,0,1], beta=[1,0,0]
                    Matrix[iMu2][iEta1] += w * get_bmoment(iMuEta, 4); //beta=[0,1,0]
                    Matrix[iMu2][iEta2] += w * get_bmoment(iMuEta, 5); //beta=[0,0,1]
                }

                if (mu1 > 0) //return iEta, iEta1, iEta2 to process the next mu1
                {
                    iEta -= eta2; //iEta -= m - eta0 +1;
                    iEta1 -= eta2;
                    iEta2 -= eta2;
                }
            }

            if (eta0 > 0) //return iMu, iMu1, iMu2 to process the next eta0
            {
                iMu -= mu2; //iMu -= n - mu0 + 1; //iMu = iMu_temp;
                iMu1 -= mu2;
                iMu2 -= mu2;
            }
        }
    }
}

void BStiff2DTri::compute_matrix(double *Fval)
{
    setFunction(Fval);
    compute_matrix();
}

void BStiff2DTri::compute_matrix(double (*f)(double, double))
{
    setFunction(f);
    compute_matrix();
}