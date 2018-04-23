#include <iostream>
using std::cin;
using std::cout;
using std::endl;
#include <cmath>

#include "Elem.h"
#include "JacobiGaussNodes.h"

BElement1D::BElement1D(int q, int n)
    : MassMat(q, n), StiffMat(q, n), LoadVec(q, n),
      ElMat(n + 1, n + 1, arma::fill::none),
      BBVector(n + 1, arma::fill::none),
      QuadVector(q, arma::fill::zeros)
{
    this->q = q;
    this->n = n;

    length = n + 1;
}

BElement1D::~BElement1D()
{
}

void BElement1D::setInterval(double a, double b)
{
    MassMat.setInterval(a, b);
    //ConvecMat.setInterval(a, b);
    StiffMat.setInterval(a, b);
    LoadVec.setInterval(a, b);
}

arma::vec BElement1D::evaluate()
{
    // initialize with 0's
    for (int i = 0; i < length; i++)
        QuadVector[i] = 0.0;

    // evaluate routine, utilizes the recurrence relation of Bernstein Polynomials
    for (int i = 0; i < length; i++)
    {
        double xi = legendre_xi(q, i);
        double s = 1 - xi;
        double r = xi / s; // recurrence relation constant
        double w = pow(s, n);
        for (int a = 0; a < length; a++)
        {
            QuadVector[i] += w * BBVector[i];
            w *= r * (n / (1. + a)); // treats the recurrence
        }
    }
    return QuadVector;
}

void BElement1D::makeSystem()
{
    MassMat.compute_matrix();
    //ConvecMat.compute_matrix();
    StiffMat.compute_matrix();
    LoadVec.compute_moments();

    for (int i = 0; i < length; i++) // length is equal to lenMass and lenStiff
        for (int j = 0; j < length; j++)
            ElMat(i, j) = MassMat.getMatrixValue(i, j) + StiffMat.getMatrixValue(i, j);
}