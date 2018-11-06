#ifndef STIFFM_H
#define STIFFM_H

#include <armadillo>
#include "Moments.h"
#define USE_DERIVATIVES

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

    // zeroes the stiffness matrix
    void zero()
    {
        BMoment1D::zero();
        Matrix.zeros();
    }

    int getLen() { return lenStiff; }

    double getMatrixValue(int i, int j)
    {
        return Matrix(i, j);
    }

    const arma::mat &getMatrix() { return Matrix; }

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

    // computes a Pascal Matrix with size lenBinomialMat
    void compute_binomials();

    // computes the normals of the triangle
    void compute_normals();

    // makes the normal vectors, unit-wise
    void normalize_normals();

    // computes all the combinations of the scalar products of the normals, storing them into N
    void compute_normals_products(arma::vec &N);

  public:
    BStiff2DTri(int q, int n);

    BStiff2DTri(int q, int n, double T[][2]);

    ~BStiff2DTri();

    // zeroes the stiffness matrix
    void zero()
    {
        BMoment2DTri::zero();
        Matrix.zeros();
    }

    int getLen() { return lenStiff; }

    double getMatrixValue(int i, int j)
    {
        return Matrix(i, j);
    }

    const arma::mat &getMatrix() { return Matrix; }

    void setFunction(const arma::vec &Fval)
    {
        BMoment2DTri::setFunction(Fval);
    }

    void setFunction(std::function<double(double, double)> f)
    {
        BMoment2DTri::setFunction(f);
    }

    void compute_matrix();

    void compute_matrix(const arma::vec &Fval);

    void compute_matrix(std::function<double(double, double)> f);
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for quadrilateral elements (2-dimensional)     *
 *****************************************************************************/
#ifdef USE_DERIVATIVES
#include "Derivatives.h"
typedef QuadD::StiffnessMatrix BStiff2DQuad;
#else
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

    const arma::mat &getMatrix() { return Matrix; }

    void setFunction(arma::vec Fval);

    void setFunction(arma::mat Fval);

    void setFunction(std::function<double(double, double)> f);

    void compute_matrix();

    void compute_matrix(arma::vec Fval);

    void compute_matrix(std::function<double(double, double)> f);
};
#endif

class BStiff3D : public BMoment3D
{
};

#endif