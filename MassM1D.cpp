#include <iostream>
using std::cin;
using std::cout;
using std::endl;

#include "MassM.h"

// maybe trade q to 2*q in the base class constructor
BMass1D::BMass1D(int q, int n) : BMoment1D(q, 2 * n)
{
    this->q = q;
    this->n = n;

    lenMass = n + 1;

    Matrix = create_matrix();
}

BMass1D::~BMass1D()
{
    delete_matrix(Matrix);
}

double **BMass1D::create_matrix()
{
    double *aux = new double[lenMass * lenMass];
    double **matrix = new double *[lenMass];
    for (int i = 0; i < lenMass; aux += lenMass, i++)
        matrix[i] = aux;
    return matrix;
}

void BMass1D::delete_matrix(double **matrix)
{
    delete matrix[0];
    delete matrix;
}

void BMass1D::compute_binomials()
{
}

void BMass1D::compute_matrix()
{
    compute_moments();
    //compute_binomials();

    for (int i = 0; i < lenMass; i++)
    {
        for (int j = 0; j < lenMass; j++)
        {
            Matrix[i][j] = (binomial(i + j, i) / binomial(2 * n, n)) * binomial(2 * n - (i + j), n - i) * get_bmoment(i);
        }
    }
}

void BMass1D::compute_matrix(double (*f)(double))
{
    setFunctionDef(f);
    compute_matrix();
}

void BMass1D::compute_matrix(double *Fval)
{
    setFunctionValue(Fval);
    compute_matrix();
}