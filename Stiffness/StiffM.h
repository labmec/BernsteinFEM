#include <armadillo>
#include "Moments.h"

#ifndef STIFFM_H
#define STIFFM_H

/*****************************************************************************
 * Bernstein Stiffness Matrix for 1-dimensional elements                     *
 *****************************************************************************/
class BStiff1D : public BMoment1D
{
    int q;                          // number of quadrature points ( recommended: 2*(n+1) )
    int n;                          // polynomial order
    int lenStiff;                   // length of the stiffness matrix
    arma::mat Matrix;               // stiffness matrix
    int lenBinomialMat;             // length of the binomial Pascal matrix
    arma::Mat<int64_t> BinomialMat; // Pascal Matrix

    // computes a Pascal Matrix with size lenBinomialMat
    void compute_binomials();

  public:
    BStiff1D(int q, int n);

    ~BStiff1D();

    int getLen() { return lenStiff; }

    double getMatrixValue(int i, int j)
    {
        return Matrix(i, j);
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

    void compute_matrix(const arma::vec &Fval);

    void compute_matrix(std::function<double(double)> f);
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for triangular elements (2-dimensional)        *
 *****************************************************************************/
class BStiff2DTri : public BMoment2DTri
{
    int q;                          // number of quadrature points ( recommended: 2*(n+1) )
    int n;                          // polynomial order
    int lenStiff;                   // length of the stiffness matrix
    arma::mat Matrix;               // stiffness matrix
    arma::mat normalMat;            // the normals of the triangle
    int lenBinomialMat;             // length of the binomial Pascal matrix
    arma::Mat<int64_t> BinomialMat; // Pascal Matrix
    BMoment2DTri *Moments;          // this is used instead of the inherited object
    // the inherited object is used to compute other parts of the matrix

    // alloc matrix linearly
    arma::mat create_matrix();

    void delete_matrix(arma::mat matrix);

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
        return Matrix(i, j);
    }

    void compute_matrix();

    void compute_matrix(arma::vec Fval);

    void compute_matrix(std::function<double(double, double)> f);
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for quadrilateral elements (2-dimensional)     *
 *****************************************************************************/
class BStiff2DQuad
{
    int q;                          // number of quadrature points ( recommended: 2*(n+1) )
    int n;                          // polynomial order
    int lenStiff;                   // length of the stiffness matrix
    arma::mat Matrix;               // stiffness matrix
    int lenBinomialMat;             // length of the binomial Pascal matrix
    BMoment2DQuad Moments1;         // moments for computation of the stiffness matrix
    BMoment2DQuad Moments2;         // moments for computation of the stiffness matrix
    arma::Mat<int64_t> BinomialMat; // Pascal Matrix

    // alloc matrix linearly
    arma::mat create_matrix();

    void delete_matrix(arma::mat matrix);

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
        return Matrix(i, j);
    }

    void setFunction(arma::vec Fval);

    void setFunction(arma::mat Fval);

    void setFunction(std::function<double(double, double)> f);

    void compute_matrix();

    void compute_matrix(arma::vec Fval);

    void compute_matrix(std::function<double(double, double)> f);
};

class BStiff3D : public BMoment3D
{
};

#endif