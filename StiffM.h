#include "Moments.h"

#ifndef STIFFM_H
#define STIFFM_H

/*****************************************************************************
 * Bernstein Stiffness Matrix for 1-dimensional elements                     *
 *****************************************************************************/
class BStiff1D : public BMoment1D
{
    int q;              // number of quadrature points ( recommended: 2*(n+1) )
    int n;              // polynomial order
    int lenStiff;       // length of the stiffness matrix
    double **Matrix;    // stiffness matrix
    int lenBinomialMat; // length of the binomial Pascal matrix
    int **BinomialMat;  // Pascal Matrix

    // alloc matrix linearly
    double **create_matrix();

    void delete_matrix(double **matrix);

    // alloc binomial matrix
    int **create_binomialMat();

    void delete_binomialMat(int **binomialMat);

    // computes a Pascal Matrix with size lenBinomialMat
    void compute_binomials();

  public:
    BStiff1D(int q, int n);

    ~BStiff1D();

    int getLen() { return lenStiff; }

    double getMatrixValue(int i, int j)
    {
        if (i < lenStiff && j < lenStiff)
            return Matrix[i][j];
        else
            return 0.0;
    }

    // computes the normalized moment
    static double grad(int k, int l)
    {
        // gradient is always 1 or -1 in the 1-dimensional case
        if (k == l)
            return 1.0;
        else
            return -1.0;
    }

    void compute_matrix();

    void compute_matrix(double *Fval);

    void compute_matrix(double (*f)(double));
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for triangular elements (2-dimensional)        *
 *****************************************************************************/
class BStiff2DTri : public BMoment2DTri
{
    int q;                  // number of quadrature points ( recommended: 2*(n+1) )
    int n;                  // polynomial order
    int lenStiff;           // length of the stiffness matrix
    double **Matrix;        // stiffness matrix
    double normalMat[3][2]; // the normals of the triangle
    int lenBinomialMat;     // length of the binomial Pascal matrix
    int **BinomialMat;      // Pascal Matrix
    BMoment2DTri *Moments;  // this is used instead of the inherited object
    // the inherited object is used to compute other parts of the matrix

    // alloc matrix linearly
    double **create_matrix();

    void delete_matrix(double **matrix);

    // alloc binomial matrix
    int **create_binomialMat();

    void delete_binomialMat(int **binomialMat);

    // computes a Pascal Matrix with size lenBinomialMat
    void compute_binomials();

    // computes the normals of the triangle
    void compute_normals();

  public:
    BStiff2DTri(int q, int n);

    BStiff2DTri(int q, int n, double T[][2]);

    ~BStiff2DTri();

    int getLen() { return lenStiff; }

    double getMatrixValue(int i, int j)
    {
        if (i < lenStiff && j < lenStiff)
            return Matrix[i][j];
        else
            return 0.0;
    }

    void compute_matrix();

    void compute_matrix(double *Fval);

    void compute_matrix(double (*f)(double, double));
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for quadrilateral elements (2-dimensional)     *
 *****************************************************************************/
class BStiff2DQuad : public BMoment2DQuad
{
    int q;              // number of quadrature points ( recommended: 2*(n+1) )
    int n;              // polynomial order
    int lenStiff;       // length of the stiffness matrix
    double **Matrix;    // stiffness matrix
    int lenBinomialMat; // length of the binomial Pascal matrix
    int **BinomialMat;  // Pascal Matrix

    // alloc matrix linearly
    double **create_matrix();

    void delete_matrix(double **matrix);

    // alloc binomial matrix
    int **create_binomialMat();

    void delete_binomialMat(int **binomialMat);

    // computes a Pascal Matrix with size lenBinomialMat
    void compute_binomials();

  public:
    BStiff2DQuad(int q, int n);

    ~BStiff2DQuad();

    int getLen() { return lenStiff; }

    double getMatrixValue(int i, int j)
    {
        if (i < lenStiff && j < lenStiff)
            return Matrix[i][j];
        else
            return 0.0;
    }

    void compute_matrix();

    void compute_matrix(double *Fval);

    void compute_matrix(double (*f)(double, double));
};

class BStiff3D : public BMoment3D
{
};

#endif