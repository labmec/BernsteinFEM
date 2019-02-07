#pragma once

#include <armadillo>
#include "Moments.h"
#include <memory>

class PascalMat;

class BMass
{
  protected:
    uint q;                          // number of quadrature points ( recommended: 2*(n+1) )
    uint n;                          // polynomial order
    uint lenMass;                    // length of the mass matrix
    arma::mat Matrix;               // mass matrix
    uint lenBinomialMat;             // length of the binomial matrix
    arma::Mat<uint64_t> BinomialMat; // Pascal Matrix

    std::shared_ptr<PascalMat> pascalMat; // Pascal matrix auxiliary obj

    void computeBinomials();

  public:
    // default constructor
    BMass(uint q, uint n);

    // copy constructor
    BMass(const BMass &cp);

    // copy assignment operator
    BMass &operator=(const BMass &cp);

    // getters

    // returns the length of one dimension of the Mass matrix
    uint length();

    // returns the polynomial order of the basis
    uint getPOrder();

    // returns a reference to the Mass Matrix
    const arma::mat &getMatrix();

    // returns the value of Matrix[i][j]
    double getMatrixValue(uint i, uint j);

  // virtual methods

    virtual void computeMatrix() = 0;
};

/*****************************************************************************
 * Bernstein Mass Matrix for 1-dimensional elements                          *
 *****************************************************************************/
class BMass1D : public BMass, public BMoment1D
{
  public:
    // default constructor
    BMass1D(uint q, uint n, const Element<Element_t::LinearEl> &el = Element<Element_t::LinearEl>());

    // copy constructor
    BMass1D(const BMass1D &cp);

    // copy assignment operator
    BMass1D &operator=(const BMass1D &cp);

    ~BMass1D();

    // compute the mass matrix
    void computeMatrix();
};

/*****************************************************************************
 * Bernstein Mass Matrix for triangular elements (2-dimensional)             *
 *****************************************************************************/
class BMass2DTri : public BMass, public BMoment2DTri
{
  public:
    // default constructor
    BMass2DTri(uint q, uint n, const Element<Element_t::TriangularEl> &el = Element<Element_t::TriangularEl>());

    // copy constructor
    BMass2DTri(const BMass2DTri &cp);

    // copy assignment operator
    BMass2DTri &operator=(const BMass2DTri &cp);

    ~BMass2DTri();

    // compute the mass matrix
    void computeMatrix();
};

/*****************************************************************************
 * Bernstein Mass Matrix for quadrilateral elements (2-dimensional)          *
 *****************************************************************************/
class BMass2DQuad : public BMass, public BMoment2DQuad
{
  public:
    // defualt constructor
    BMass2DQuad(uint q, uint n, Element<Element_t::QuadrilateralEl> element = Element<Element_t::QuadrilateralEl>());

    // copy constructor
    BMass2DQuad(const BMass2DQuad &cp);

    //copy assignment operator
    BMass2DQuad &operator=(const BMass2DQuad &cp);

    ~BMass2DQuad();

    // compute the mass matrix
    void computeMatrix();
};

class BMass3DCube : public BMoment3DCube
{
};