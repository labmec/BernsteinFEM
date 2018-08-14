#include "Derivatives.h"
#include "JacobiGaussNodes.h"

#include <iostream>

#ifndef LEN
#define LEN(N) ((N + 1) * (N + 1))
#endif

using namespace QuadD;
using namespace arma;

void QuadDerivative::compute_binomials(Mat<int64_t> &BinomialMat, int lenBinom)
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

QuadDerivative::QuadDerivative(int q, int n)
    : Matrix(LEN(n), LEN(n), fill::zeros),
      BinomialMat(n + 1, n + 1, fill::zeros),
      Fval(q * q, fill::none)

{
    this->q = q;
    this->n = n;

    len = LEN(n);
    lenBinom = n + 1;

    compute_binomials(BinomialMat, lenBinom);
}

void QuadDerivative::setFunction(const vec &Fval)
{
    this->Fval = Fval;
}

arma::mat QuadDerivative::getIntegrationPoints()
{
    int q = QuadDerivative::q;
    arma::mat points(q * q, 2);
    arma::vec X(2, arma::fill::none);
    BMoment2DQuad nodalAux(1, 1);
    nodalAux.setQuadrilateral(this->vertices);
    double dX;

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            double xi = (legendre_xi(q, i) + 1.0) * 0.5;
            double eta = (legendre_xi(q, j) + 1.0) * 0.5;
            nodalAux.nodalShape(X, dX, xi, eta);
            points(i * q + j, 0) = X[0];
            points(i * q + j, 1) = X[1];
        }
    }
    
    return points;
}

/****************************************
 ********** Defining dXi_dXi ************
 ****************************************/

dXi_dXi::dXi_dXi(int q, int n)
    : QuadDerivative(q, n) { }

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
                    double w = wXi * wEta * Fval(iXi * q + iEta);
                    moments(aXi, aEta) += w;

                    wEta *= rEta * (nEta - aEta) / (1.0 + aEta);
                }
            }
            wXi *= rXi * (nXi - aXi) / (1.0 + aXi);
        }
    } // O(q^2 * n^2)


    // test if the moments were calculated correctly
    // std::ofstream file;
    // file.open("results/xi_xi_moments.txt");
    // moments.print(file);
    // file.close();

    // then we arrange the terms

    double Const = n * n * (1.0 / (BinomialMat(n, n) * BinomialMat(n - 1, n - 1)));

    for (int a1 = 0; a1 < n; a1++)
    {
        for (int b1 = 0; b1 < n; b1++)
        {
            double w1 = Const * BinomialMat(a1, b1) * BinomialMat(n - a1 - 1, n - b1 - 1);
            for (int a2 = 0; a2 <= n; a2++)
            {
                for (int b2 = 0; b2 <= n; b2++)
                {
                    double w2 = w1 * BinomialMat(a2, b2) * BinomialMat(n - a2, n - b2);
                    double mom = w2 * moments(a1 + b1, a2 + b2);

                    int i = BMoment2DQuad::position(a1, a2, n);
                    int j = BMoment2DQuad::position(b1, b2, n);
                    int I = BMoment2DQuad::position(a1 + 1, a2, n);
                    int J = BMoment2DQuad::position(b1 + 1, b2, n);

                    Matrix(i, j) += mom;

                    Matrix(I, J) += mom;

                    Matrix(I, j) -= mom;

                    Matrix(i, J) -= mom;
                }
            }
        }
    } // O(n^4)
}

/******************************************
 ********** Defining dEta_dEta ************
 ******************************************/

dEta_dEta::dEta_dEta(int q, int n)
    : QuadDerivative(q, n) { }

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
        double sXi = 1.0 - xi;
        double rXi = xi / sXi;
        double omegaXi = legendre_w(q, iXi) * 0.5;

        double wXi = omegaXi * pow(sXi, nXi);
        for (int aXi = 0; aXi <= nXi; aXi++)
        {
            for (int iEta = 0; iEta < q; iEta++)
            {
                double eta = (legendre_xi(q, iEta) + 1.0) * 0.5;
                double sEta = 1.0 - eta;
                double rEta = eta / sEta;
                double omegaEta = legendre_w(q, iEta) * 0.5;

                double wEta = omegaEta * pow(sEta, nEta);
                for (int aEta = 0; aEta <= nEta; aEta++)
                {
                    moments(aXi, aEta) += wXi * wEta * Fval(iXi * q + iEta);

                    wEta *= rEta * ((nEta - aEta) / (1.0 + aEta));
                }
            }
            wXi *= rXi * ((nXi - aXi) / (1.0 + aXi));
        }
    } // O(q^2 * n^2)

    // test if the moments were calculated correctly
    // std::ofstream file;
    // file.open("results/eta_eta_moments.txt");
    // moments.print(file);
    // file.close();

    // then we arrange the terms

    double Const = n * n * (1.0 / (BinomialMat(n, n) * BinomialMat(n - 1, n - 1)));

    for (int a1 = 0; a1 <= n; a1++)
    {
        for (int b1 = 0; b1 <= n; b1++)
        {
            double w1 = Const * BinomialMat(a1, b1) * BinomialMat(n - a1, n - b1);
            for (int a2 = 0; a2 < n; a2++)
            {
                for (int b2 = 0; b2 < n; b2++)
                {
                    double w2 = w1 * BinomialMat(a2, b2) * BinomialMat(n - a2 - 1, n - b2 - 1);
                    double mom = w2 * moments(a1 + b1, a2 + b2);

                    int i = BMoment2DQuad::position(a1, a2, n);
                    int j = BMoment2DQuad::position(b1, b2, n);
                    int I = BMoment2DQuad::position(a1, a2 + 1, n);
                    int J = BMoment2DQuad::position(b1, b2 + 1, n);

                    Matrix(i, j) += mom;

                    Matrix(I, J) += mom;

                    Matrix(I, j) -= mom;

                    Matrix(i, J) -= mom;
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
      QuadDerivative(q, n) { }

void dXi_dEta::compute_matrix()
{
    compute_moments();

    // test if the moments were calculated correctly
    // std::ofstream file;
    // file.open("results/xi_eta_moments.txt");
    // BMoment2DQuad::get_bmoment().print(file);
    // file.close();

    Matrix.zeros();
    int n = QuadDerivative::n; // takes the n from QuadDerivative, the BMoment2DQuad n is equal to this one times 2 minus 1

    // then arrange the terms
    double Const = n * n * (1.0 / (BinomialMat(n - 1, n) * BinomialMat(n, n - 1))); // don't worry the BinomialMat is symmetrical

    for (int a1 = 0; a1 < n; a1++)
    {
        for (int b1 = 0; b1 <= n; b1++)
        {
            double w1 = Const * BinomialMat(a1, b1) * BinomialMat(n - a1 - 1, n - b1);
            for (int a2 = 0; a2 <= n; a2++)
            {
                for (int b2 = 0; b2 < n; b2++)
                {
                    double w2 = w1 * BinomialMat(a2, b2) * BinomialMat(n - a2, n - b2 - 1);
                    double mom = w2 * get_bmoment((a1 + b1) * (n + n) + a2 + b2);

                    int i = position(a1, a2, n);
                    int j = position(b1, b2, n);
                    int I = position(a1 + 1, a2, n); // a1 + 1
                    int J = position(b1, b2 + 1, n); // b2 + 1

                    Matrix(i, j) += mom;
                    Matrix(I, J) += mom;

                    Matrix(I, j) -= mom;
                    Matrix(i, J) -= mom;
                }
            }
        }
    }
} // takes the time to compute moments with parameters q and 2 * n - 2 (q^2 * (2 * n - 2)^2) plus n^4

/****************************************
 ********** Stiffness Matrix ************
 ****************************************/

StiffnessMatrix::StiffnessMatrix(int q, int n)
    : QuadDerivative(q, n), Xi_Xi(q, n), Xi_Eta(q, n), Eta_Eta(q, n) {}

void StiffnessMatrix::setFunction(const vec &Fval)
{
    // auxiliary object to compute nodal shape function
    BMoment2DQuad nodalAux(1, 1);
    nodalAux.setQuadrilateral(this->vertices);
    mat jac(2, 2, fill::none);
    vec X(2, fill::none);

    vec Fval1 = Fval; // Xi_Xi function values
    vec Fval2 = Fval; // Eta_Eta function values
    vec Fval3 = Fval; // Xi_eta function values

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            double dX;
            double xi = (legendre_xi(q, i) + 1.0) * 0.5;
            double eta = (legendre_xi(q, j) + 1.0) * 0.5;
            nodalAux.nodalShape(X, jac, dX, xi, eta);
            mat jac_inv = jac.i();

            Fval1(i * q + j) *= (jac_inv(0, 0) * jac_inv(0, 0) + jac_inv(0, 1) * jac_inv(0, 1)) * dX;
            Fval2(i * q + j) *= (jac_inv(1, 0) * jac_inv(1, 0) + jac_inv(1, 1) * jac_inv(1, 1)) * dX;
            Fval3(i * q + j) *= (jac_inv(0, 0) * jac_inv(1, 0) + jac_inv(0, 1) * jac_inv(1, 1)) * dX;
        }
    }

    Xi_Xi.setFunction(Fval1);
    Eta_Eta.setFunction(Fval2);
    Xi_Eta.setFunction(Fval3);
}

void StiffnessMatrix::setQuadrilateral(const mat &vertices)
{
    this->vertices = vertices;
}

// computes the stiffness matrix from the coefficients of partial derivatives in the master element
void StiffnessMatrix::compute_matrix()
{
    Xi_Xi.compute_matrix();
    Xi_Eta.compute_matrix();
    Eta_Eta.compute_matrix();

    mat xi_xi = Xi_Xi.getMatrix();
    mat xi_eta = Xi_Eta.getMatrix();
    mat eta_eta = Eta_Eta.getMatrix();

    std::ofstream f1, f2, f3;
    f1.open("results/xi_xi.txt");
    f2.open("results/xi_eta.txt");
    f3.open("results/eta_eta.txt");
    xi_xi.print(f1);
    xi_eta.print(f2);
    eta_eta.print(f3);
    f1.close();
    f2.close();
    f3.close();

    for (int i = 0; i < LEN(n); i++)
    {
        for (int j = 0; j < LEN(n); j++)
        {
            Matrix(i, j) = xi_xi(i, j) + eta_eta(i, j) + xi_eta(i, j) + xi_eta(j, i);
        }
    }
}