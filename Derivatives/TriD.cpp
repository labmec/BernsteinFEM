#include "Derivatives.h"
#include "JacobiGaussNodes.h"

#ifndef LEN
#define LEN(N) ((MAX(n, q - 1) + 1) * (MAX(n, q - 1) + 1))
#endif

using namespace TriD;
using namespace arma;

TriangleDerivative::TriangleDerivative(int q, int n)
    : Matrix(LEN(n), LEN(n), fill::zeros),
      BinomialMat(n + 1, n + 1, fill::zeros),
      Fval(q * q, fill::none)
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

double TriangleDerivative::area()
{
    double x1 = vertices(0, 0);
    double y1 = vertices(0, 1);
    double x2 = vertices(1, 0);
    double y2 = vertices(1, 1);
    double x3 = vertices(2, 0);
    double y3 = vertices(2, 1);

    return abs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2.;
}

void TriangleDerivative::setFunction(const mat &Fval)
{
    this->Fval = Fval;
}

void TriangleDerivative::compute_binomials(Mat<int64_t> &BinomialMat, int lenBinom)
{
    for (int i = 0; i < lenBinom; i++)
    {
        BinomialMat.at(i, 0) = 1;
        BinomialMat.at(0, i) = 1;
    }

    for (int k = 1; k < lenBinom; k++)
    {
        for (int l = 1; l < lenBinom; l++)
        {
            BinomialMat.at(k, l) += BinomialMat.at(k, l - 1) + BinomialMat.at(k - 1, l);
        }
    }
}

void TriangleDerivative::print_mathematica(std::ostream &stream)
{
    stream.flags(std::ios_base::fixed);
    stream.precision(16);
    stream << "{";
    for (int i = 0; i < (n + 1) * (n + 1); i++)
    {
        stream << "{";
        for (int j = 0; j < (n + 1) * (n + 1); j++)
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

void TriangleDerivative::print(std::ostream &stream)
{
    for (int i = 0; i < (n + 1) * (n + 1); i++)
    {
        for (int j = 0; j < (n + 1) * (n + 1); j++)
            stream << std::scientific << Matrix(i, j) << "\t";
        stream << "\n";
    }
}


const double &TriangleDerivative::operator()(int i, int j)
{
    try
    {
        return Matrix(i, j);
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "i = " << i << std::endl
                  << "j = " << j << std::endl
                  << "Max: " << len << std::endl;
        abort();
    }
}

/****************************************
 ********** Defining dXi_dXi ************
 ****************************************/

dXi_dXi::dXi_dXi(int q, int n)
    : TriangleDerivative(q, n) {}

void dXi_dXi::compute_matrix()
{
    int N = 2 * n - 2;

    vec mu_0((N + 1) * (N + 1), fill::zeros);
    vec mu_1((N + 1) * (N + 3), fill::zeros);
    vec mu_2((N + 1) * (N + 2), fill::zeros);

    vec mu_0_inter((N + 1) * q, fill::zeros);
    vec mu_1_inter((N + 1) * q, fill::zeros);
    vec mu_2_inter((N + 1) * q, fill::zeros);

    // calculate mu_0, mu_1, mu_2
    // no scaling constant, for they are already in the jacobian of the Stiffness matrix

    // convert first index for all moments
    for (int i = 0; i < q; i++)
    {
        double xi = (1.0 + legendre_xi(q, i)) * 0.5;
        double wgt = legendre_w(q, i) * 0.5;

        double s = 1 - xi;
        double r = xi / s;

        double B = wgt * pow(s, N);
        for (int a1 = 0; a1 <= N; a1++)
        {
            for (int j = 0; j < q; j++)
            {
                double aux = B * Fval(i * q + j);
                mu_0_inter(a1 * q + j) += aux;
                mu_1_inter(a1 * q + j) += aux;
                mu_2_inter(a1 * q + j) += aux;
            }
            B = B * r * (N - a1) / (1 + a1);
        }
    }

    mu_0_inter.print("mu_0_inter:");
    mu_1_inter.print("mu_1_inter:");
    mu_2_inter.print("mu_2_inter:");

    // convert second index for all moments
    for (int i = 0; i < q; i++)
    {
        double xi = (1.0 + legendre_xi(q, i)) * 0.5;
        double wgt = legendre_w(q, i) * 0.5;

        double s = 1 - xi;
        double r = xi / s;

        for (int a1 = 0; a1 <= N; a1++)
        {
            double B0 = wgt * pow(s, N - a1);
            double B2 = B0 * s; // s^(N-a1+1)
            double B1 = B2 * s; // s^(N-a1+2)
            double m0 = mu_0_inter(a1 * q + i);
            double m1 = mu_1_inter(a1 * q + i);
            double m2 = mu_2_inter(a1 * q + i);
            for (int a2 = 0; a2 <= N - a1; a2++)
            {
                mu_0(a1 * (N + 1) + a2) += B0 * m0;
                mu_1(a1 * (N + 3) + a2) += B1 * m1;
                mu_2(a1 * (N + 2) + a2) += B2 * m2;

                B0 = B0 * r * (N - a1 - a2) / (1 + a2);
                B1 = B1 * r * (N - a1 - a2 + 2) / (1 + a2);
                B2 = B2 * r * (N - a1 - a2 + 1) / (1 + a2);
            }
            mu_1(a1 * (N + 3) + N - a1 + 1) += B1 * m1;
            B1 = B1 * r / (1 + (N - a1 + 1));
            mu_1(a1 * (N + 3) + N - a1 + 2) += B1 * m1;
            mu_2(a1 * (N + 2) + N - a1 + 1) += B2 * m2;
        }
    }

    mu_0.print("mu_0:");
    mu_1.print("mu_1:");
    mu_2.print("mu_2:");
    double Const = n * n * (1.0 / BinomialMat(n - 1, n - 1));

    // compute matrix using mu_0, mu_1 and mu_2 in O(n^4)
    for (int a1 = 0; a1 < n - 1; a1++)
    {
        for (int b1 = 0; b1 < n - 1; b1++)
        {
            double w1 = Const * BinomialMat(a1, b1);
            for (int a2 = 0; a2 <= n - a1; a2++)
            {
                for (int b2 = 0; b2 <= n - b1; b2++)
                {
                    double w2 = w1 * BinomialMat(a2, b2) * BinomialMat(n - a1 - a2, n - b1 - b2);

                    int i = BMoment2DTri::position(a1, a2, n);
                    int j = BMoment2DTri::position(b1, b2, n);
                    int I = BMoment2DTri::position(a1 + 1, a2, n);
                    int J = BMoment2DTri::position(b1 + 1, b2, n);

                    Matrix(i, j) += w2 * mu_1((a1 + b1) * (n + 3) + a2 + b2);
                    Matrix(I, J) += w2 * mu_0((a1 + b1) * (n + 1) + a2 + b2);
                    Matrix(i, J) -= w2 * mu_2((a1 + b1) * (n + 2) + a2 + b2);
                    Matrix(I, j) -= w2 * mu_2((a1 + b1) * (n + 2) + a2 + b2);
                }
            }
        }
    }
}

/*****************************************
 ********** Defining dEta_dEta ***********
 *****************************************/

dEta_dEta::dEta_dEta(int q, int n)
    : TriangleDerivative(q, n) {}

void dEta_dEta::compute_matrix()
{
}

/****************************************
 ********** Defining dXi_dEta ***********
 ****************************************/

dXi_dEta::dXi_dEta(int q, int n)
    : TriangleDerivative(q, n) {}

void dXi_dEta::compute_matrix()
{
}