#include "Derivatives.h"
#include "JacobiGaussNodes.h"

#include <iostream>

#ifndef LEN
#define LEN(N) ((N + 1) * (N + 1))
#endif

using namespace QuadD;
using namespace arma;

void QuadDerivative::compute_binomials(Mat<int64_t> &BinomialMat, uint lenBinom)
{
    for (uint i = 0; i < lenBinom; i++)
    {
        BinomialMat.at(i, 0) = 1;
        BinomialMat.at(0, i) = 1;
    }

    for (uint k = 1; k < lenBinom; k++)
    {
        for (uint l = 1; l < lenBinom; l++)
        {
            BinomialMat.at(k, l) += BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

QuadDerivative::QuadDerivative(uint q, uint n)
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
    arma::mat points(q * q, 2);

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            points.at(i * q + j, 0) = (legendre_xi(q, i) + 1.0) * 0.5;
            points.at(i * q + j, 1) = (legendre_xi(q, j) + 1.0) * 0.5;
        }
    }

    return points;
}

void QuadDerivative::print_mathematica(std::ostream &stream)
{
    stream.flags(std::ios_base::fixed);
    stream.precision(16);
    stream << "{";
    for (uint i = 0; i < (n + 1) * (n + 1); i++)
    {
        stream << "{";
        for (uint j = 0; j < (n + 1) * (n + 1); j++)
            if (j < (n + 1) * (n + 1) - 1)
                stream << Matrix(i, j) << ", ";
            else
                stream << Matrix(i, j);
        if (i < (n + 1) * (n + 1) - 1)
            stream << "},\n";
        else
            stream << "}";
    }
    stream << "}" << endl;
}

void QuadDerivative::print(std::ostream &stream)
{
    for (uint i = 0; i < (n + 1) * (n + 1); i++)
    {
        for (uint j = 0; j < (n + 1) * (n + 1); j++)
            stream << std::scientific << Matrix(i, j) << "\t";
        stream << "\n";
    }
}

/****************************************
 ********** Defining dXi_dXi ************
 ****************************************/

dXi_dXi::dXi_dXi(uint q, uint n)
    : BMoment2DQuad(q, 2 * (n - 1), 2 * n),
      QuadDerivative(q, n) {}

void dXi_dXi::compute_matrix()
{
    BMoment2DQuad::computeMoments();

    // test if the moments were calculated correctly
    // std::ofstream file;
    // file.open("results/xi_xi_moments.txt");
    // moments.print(file);
    // file.close();

    // then we arrange the terms

    uint n = QuadDerivative::n;

    double Const = n * n * (1.0 / (BinomialMat.at(n, n) * BinomialMat.at(n - 1, n - 1)));

    for (uint a1 = 0; a1 < n; a1++)
    {
        for (uint b1 = 0; b1 < n; b1++)
        {
            double w1 = Const * BinomialMat.at(a1, b1) * BinomialMat.at(n - a1 - 1, n - b1 - 1);
            for (uint a2 = 0; a2 <= n; a2++)
            {
                for (uint b2 = 0; b2 <= n; b2++)
                {
                    double w2 = w1 * BinomialMat.at(a2, b2) * BinomialMat.at(n - a2, n - b2);
                    double mom = w2 * getBMoment(a1 + b1, a2 + b2, 0);

                    uint i = element.position({a1, a2}, n);
                    uint j = element.position({b1, b2}, n);
                    uint I = element.position({a1 + 1, a2}, n);
                    uint J = element.position({b1 + 1, b2}, n);

                    Matrix.at(i, j) += mom;

                    Matrix.at(I, J) += mom;

                    Matrix.at(I, j) -= mom;

                    Matrix.at(i, J) -= mom;
                }
            }
        }
    } // O(n^4)
}

/******************************************
 ********** Defining dEta_dEta ************
 ******************************************/

dEta_dEta::dEta_dEta(uint q, uint n)
    : BMoment2DQuad(q, 2 * n, 2 * (n - 1)),
      QuadDerivative(q, n) {}

void dEta_dEta::compute_matrix()
{
    BMoment2DQuad::computeMoments();
    // dEta_dEta is basically the same as dXi_dXi, we could do it by transposing the FVal matrix
    // and then transposing the resulting matrix

    uint n = QuadDerivative::n;

    // test if the moments were calculated correctly
    // std::ofstream file;
    // file.open("results/eta_eta_moments.txt");
    // moments.print(file);
    // file.close();

    // then we arrange the terms

    double Const = n * n * (1.0 / (BinomialMat.at(n, n) * BinomialMat.at(n - 1, n - 1)));

    for (uint a1 = 0; a1 <= n; a1++)
    {
        for (uint b1 = 0; b1 <= n; b1++)
        {
            double w1 = Const * BinomialMat.at(a1, b1) * BinomialMat.at(n - a1, n - b1);
            for (uint a2 = 0; a2 < n; a2++)
            {
                for (uint b2 = 0; b2 < n; b2++)
                {
                    double w2 = w1 * BinomialMat.at(a2, b2) * BinomialMat.at(n - a2 - 1, n - b2 - 1);
                    double mom = w2 * getBMoment(a1 + b1, a2 + b2, 0);

                    uint i = element.position({a1, a2}, n);
                    uint j = element.position({b1, b2}, n);
                    uint I = element.position({a1, a2 + 1}, n);
                    uint J = element.position({b1, b2 + 1}, n);

                    Matrix.at(i, j) += mom;

                    Matrix.at(I, J) += mom;

                    Matrix.at(I, j) -= mom;

                    Matrix.at(i, J) -= mom;
                }
            }
        }
    } // O(n^4)
}

/****************************************
 ********** Defining dXi_dEta ************
 ****************************************/

dXi_dEta::dXi_dEta(uint q, uint n)
    : BMoment2DQuad(q, 2 * n - 1),
      QuadDerivative(q, n) {}

void dXi_dEta::compute_matrix()
{
    BMoment2DQuad::computeMoments();

    // test if the moments were calculated correctly
    // std::ofstream file;
    // file.open("results/xi_eta_moments.txt");
    // BMoment2DQuad::get_bmoment().print(file);
    // file.close();

    uint n = QuadDerivative::n; // takes the n from QuadDerivative, the BMoment2DQuad n is equal to this one times 2 minus 1

    // then arrange the terms
    double Const = n * n * (1.0 / (BinomialMat.at(n - 1, n) * BinomialMat.at(n, n - 1))); // don't worry the BinomialMat is symmetrical

    for (uint a1 = 0; a1 < n; a1++)
    {
        for (uint b1 = 0; b1 <= n; b1++)
        {
            double w1 = Const * BinomialMat.at(a1, b1) * BinomialMat.at(n - a1 - 1, n - b1);
            for (uint a2 = 0; a2 <= n; a2++)
            {
                for (uint b2 = 0; b2 < n; b2++)
                {
                    double w2 = w1 * BinomialMat.at(a2, b2) * BinomialMat.at(n - a2, n - b2 - 1);
                    double mom = w2 * getBMoment((a1 + b1) * (n + n) + a2 + b2);

                    int i = element.position({a1, a2}, n);
                    int j = element.position({b1, b2}, n);
                    int I = element.position({a1 + 1, a2}, n); // a1 + 1
                    int J = element.position({b1, b2 + 1}, n); // b2 + 1

                    Matrix.at(i, j) += mom;
                    Matrix.at(I, J) += mom;

                    Matrix.at(I, j) -= mom;
                    Matrix.at(i, J) -= mom;
                }
            }
        }
    }
} // takes the time to compute moments with parameters q and 2 * n - 2 (q^2 * (2 * n - 2)^2) plus n^4

/****************************************
 ********** Stiffness Matrix ************
 ****************************************/

StiffnessMatrix::StiffnessMatrix(uint q, uint n, const Element<Element_t::QuadrilateralEl> &el)
    : QuadDerivative(q, n), Xi_Xi(q, n), Xi_Eta(q, n), Eta_Eta(q, n), element(el) {}

void StiffnessMatrix::setFunction(const vec &Fval)
{
    // auxiliary object to compute nodal shape function
    mat jac(2, 2, fill::none);
    vec X(2, fill::none);

    vec Fval1 = Fval; // Xi_Xi function values
    vec Fval2 = Fval; // Eta_Eta function values
    vec Fval3 = Fval; // Xi_eta function values

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            double xi = (legendre_xi(q, i) + 1.0) * 0.5;
            double eta = (legendre_xi(q, j) + 1.0) * 0.5;
            element.mapToElement({xi, eta}, jac);
            double dX = det(jac);
            mat jac_inv = inv(jac);

            Fval1.at(i * q + j) *= (jac_inv.at(0, 0) * jac_inv.at(0, 0) + jac_inv.at(0, 1) * jac_inv.at(0, 1)) * dX;
            Fval2.at(i * q + j) *= (jac_inv.at(1, 0) * jac_inv.at(1, 0) + jac_inv.at(1, 1) * jac_inv.at(1, 1)) * dX;
            Fval3.at(i * q + j) *= (jac_inv.at(0, 0) * jac_inv.at(1, 0) + jac_inv.at(0, 1) * jac_inv.at(1, 1)) * dX;
        }
    }

    Xi_Xi.setFunction(Fval1);
    Eta_Eta.setFunction(Fval2);
    Xi_Eta.setFunction(Fval3);
}

void StiffnessMatrix::setQuadrilateral(const mat &vertices)
{
    element = Element<Element_t::QuadrilateralEl>(vertices);
}

arma::mat StiffnessMatrix::getIntegrationPoints()
{
    uint q = QuadDerivative::q;
    arma::mat points(q * q, 2);

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            double xi = (legendre_xi(q, i) + 1.0) * 0.5;
            double eta = (legendre_xi(q, j) + 1.0) * 0.5;
            auto X = element.mapToElement({xi, eta});
            points.at(i * q + j, 0) = X[0];
            points.at(i * q + j, 1) = X[1];
        }
    }

    return points;
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

    // std::ofstream f1, f2, f3;
    // f1.open("results/xi_xi.txt");
    // f2.open("results/xi_eta.txt");
    // f3.open("results/eta_eta.txt");
    // xi_xi.print(f1);
    // xi_eta.print(f2);
    // eta_eta.print(f3);
    // f1.close();
    // f2.close();
    // f3.close();

    for (uint i = 0; i < LEN(n); i++)
    {
        for (uint j = 0; j < LEN(n); j++)
        {
            Matrix.at(i, j) = xi_xi.at(i, j) + eta_eta.at(i, j) + xi_eta.at(i, j) + xi_eta.at(j, i);
        }
    }
}