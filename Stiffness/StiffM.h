#pragma once

#include <armadillo>
#include "Moments.h"
#define USE_DERIVATIVES

class BStiff
{
  protected:
    int q;                          // number of quadrature points ( recommended: 2*(n+1) )
    int n;                          // polynomial order
    int lenStiff;                   // length of the stiffness matrix
    arma::mat Matrix;               // stiffness matrix
    int lenBinomialMat;             // length of the binomial Pascal matrix
    arma::Mat<int64_t> BinomialMat; // Pascal Matrix

    // computes a Pascal Matrix with size lenBinomialMat
    void computeBinomials();

  public:
    // default constructor
    BStiff(int q, int n);

    // copy constructor
    BStiff(const BStiff &cp);

    // copy assignment operator (unnecessary)
    BStiff &operator=(const BStiff &cp);

  // getters

    // returns the length of the matrix
    int length();

    // returns the polynomial order of the basis
    int getPOrder();

    // returns the value of Matrix[i][j]
    double getMatrixValue(int i, int j);

    // returns the matrix
    const arma::mat &getMatrix();

  // other methods

    // assign 0 to all elements of the matrix
    virtual void zero();

  // abstract methods
    
    // compute the stiffness matrix
    virtual void computeMatrix() = 0;
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for 1-dimensional elements                     *
 *****************************************************************************/
class BStiff1D : public BStiff, public BMoment1D
{
    // computes the normalized moment, 
    static double grad(int k, int l)
    {
        // gradient is always 1 or -1 in the 1-dimensional case
        if (k == l)
            return 1.0;
        else
            return -1.0;
    }
  public:
    // default constructor
    BStiff1D(int q, int n, const Element<Element_t::LinearEl> &element = Element<Element_t::LinearEl>());

    // copy constructor
    BStiff1D(const BStiff1D &cp);

    // copy assignment operator
    BStiff1D &operator=(const BStiff1D &cp);

    ~BStiff1D();

    // zeroes the stiffness matrix
    inline void zero()
    {
        BStiff::zero();
        BMoment1D::zero();
    }

    void computeMatrix();
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for triangular elements (2-dimensional)        *
 *****************************************************************************/
class BStiff2DTri : public BStiff, public BMoment2DTri
{
    arma::mat normalMat; // normal vectors to each triangle side

    // computes the normals of the triangle
    void compute_normals();

    // makes the normal vectors, unit-wise
    void normalize_normals();

    // computes all the combinations of the scalar products of the normals, storing them into N
    void compute_normals_products(arma::vec &N);

  public:
    // default constructor
    BStiff2DTri(int q, int n, const Element<Element_t::TriangularEl> &element = Element<Element_t::TriangularEl>());

    // copy constructor
    BStiff2DTri(const BStiff2DTri &cp);

    // copy assignment operator
    BStiff2DTri &operator=(const BStiff2DTri &cp);

    ~BStiff2DTri();

    // zeroes the stiffness matrix
    inline void zero()
    {
        BStiff::zero();
        BMoment2DTri::zero();
    }

    void computeMatrix();
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
  public:
    // default constructor
    BStiff2DQuad(int q, int n, const Element<Element_t::QuadrilateralEl> &element = Element<Element_t::QuadrilateralEl>());

    // copy constructor
    BStiff2DQuad(const BStiff2DQuad &cp);

    // copy assignment operator
    BSitff2DQuad &operator=(const BStiff2DQuad &cp);

    ~BStiff2DQuad();

    void computeMatrix();
};
#endif

class BStiff3DCube : public BMomentCube3D
{
};