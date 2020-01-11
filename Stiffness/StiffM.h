#pragma once

#include "Moments.h"
#define USE_DERIVATIVES

class BStiff
{
  protected:
    uint32_t q;                          // number of quadrature points ( recommended: 2*(n+1) )
    uint32_t n;                          // polynomial order
    uint32_t lenStiff;                   // length of the stiffness matrix
    TPZFMatrix<REAL> Matrix;               // stiffness matrix
    uint32_t lenBinomialMat;             // length of the binomial Pascal matrix
    TPZFMatrix<int64_t> BinomialMat; // Pascal Matrix

    // computes a Pascal Matrix with size lenBinomialMat
    void computeBinomials();

  public:
    // default constructor
    BStiff(uint32_t q, uint32_t n);

    // copy constructor
    BStiff(const BStiff &cp);

    // copy assignment operator (unnecessary)
    BStiff &operator=(const BStiff &cp);

  // getters

    // returns the length of the matrix
    uint32_t length();

    // returns the polynomial order of the basis
    uint32_t getPOrder();

    // returns the value of Matrix[i][j]
    REAL getMatrixValue(uint32_t i, uint32_t j);

    // returns the matrix
    const TPZFMatrix<REAL> &getMatrix();

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
  public:
    // default constructor
    BStiff1D(uint32_t q, uint32_t n, const Element<Element_t::LinearEl> &element = Element<Element_t::LinearEl>());

    // copy constructor
    BStiff1D(const BStiff1D &cp);

    // copy assignment operator
    BStiff1D &operator=(const BStiff1D &cp);

    virtual ~BStiff1D() = default;

    // zeroes the stiffness matrix
    inline void zero() override
    {
        BStiff::zero();
        BMoment1D::zero();
    }

    void computeMatrix() override;
};

/*****************************************************************************
 * Bernstein Stiffness Matrix for triangular elements (2-dimensional)        *
 *****************************************************************************/
class BStiff2DTri : public BStiff, public BMoment2DTri
{
    TPZFMatrix<REAL> normalMat; // normal vectors to each triangle side

    // computes the normals of the triangle
    void compute_normals();

    // makes the normal vectors, unit-wise
    void normalize_normals();

    // computes all the combinations of the scalar products of the normals, storing them into N
    void compute_normals_products(TPZVec<REAL> &N);

  public:
    // default constructor
    BStiff2DTri(uint32_t q, uint32_t n, const Element<Element_t::TriangularEl> &element = Element<Element_t::TriangularEl>());

    // copy constructor
    BStiff2DTri(const BStiff2DTri &cp);

    // copy assignment operator
    BStiff2DTri &operator=(const BStiff2DTri &cp);

    virtual ~BStiff2DTri() = default;

    // zeroes the stiffness matrix
    inline void zero() override
    {
        BStiff::zero();
        BMoment2DTri::zero();
    }

    void computeMatrix() override;
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
    BStiff2DQuad(uint32_t q, uint32_t n, const Element<Element_t::QuadrilateralEl> &element = Element<Element_t::QuadrilateralEl>());

    // copy constructor
    BStiff2DQuad(const BStiff2DQuad &cp);

    // copy assignment operator
    BSitff2DQuad &operator=(const BStiff2DQuad &cp);

    ~BStiff2DQuad();

    void computeMatrix();
};
#endif

class BStiff3DCube : public BMoment3DCube
{
};

class BStiff3DTetra : public BMoment3DTetra
{

};

// BStiff copy factory
template<Element_t EL>
BStiff* copyBStiff(BStiff const& copy);