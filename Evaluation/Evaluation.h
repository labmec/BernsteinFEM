#pragma once

#include "Elements.h"
#include <armadillo>

// abstract class
template <Element_t EL>
class BEval
{
  protected:
    int q;               // number of integration points
    int n;               // polynomial order
    arma::vec BBVec;     // Vector holding coefficients
    uint bbvec_len;      // length BBVec should have
    arma::vec eval;      // vector holding the polynomial evaluated in the integration points of the element
    Element<EL> element; // element of computation

  public:
    BEval(int q, int n, Element<EL> const &element = Element<EL>());

    BEval(int q, int n, arma::vec const &coeffVec, Element<EL> const &element = Element<EL>());

    // getters

    int getNumIntegrationPoints();

    int getPOrder();

    arma::vec &getCoefficientsVec();

    arma::vec &getEval();

    Element<EL> &getElement();

    // setters

    void setCoefficientsVec(arma::vec const &vec);

    // evaluates the resulting polynomial at the integration points
    // and returns the resulting vector
    virtual arma::vec &computeEvaluation() = 0;
};

class BEval1D : public BEval<Element_t::LinearEl>
{
  public:
    BEval1D(int q, int n, Element<Element_t::LinearEl> const &element = Element<Element_t::LinearEl>());

    BEval1D(int q, int n, arma::vec const &coeffVec, Element<Element_t::LinearEl> const &element = Element<Element_t::LinearEl>());

    arma::vec &computeEvaluation();
};

class BEval2DTri : public BEval<Element_t::TriangularEl>
{
    arma::vec eval_inter;

  public:
    BEval2DTri(int q, int n, Element<Element_t::TriangularEl> const &element = Element<Element_t::TriangularEl>());

    BEval2DTri(int q, int n, arma::vec const &coeffVec, Element<Element_t::TriangularEl> const &element = Element<Element_t::TriangularEl>());

    arma::vec &computeEvaluation();
};

class BEval2DQuad : public BEval<Element_t::QuadrilateralEl>
{
    arma::vec eval_inter;

  public:
    BEval2DQuad(int q, int n, Element<Element_t::QuadrilateralEl> const &element = Element<Element_t::QuadrilateralEl>());

    BEval2DQuad(int q, int n, arma::vec const &coeffVec, Element<Element_t::QuadrilateralEl> const &element = Element<Element_t::QuadrilateralEl>());

    arma::vec &computeEvaluation();
};

class BEval3DCube : public BEval<Element_t::CubeEl>
{
  public:
    BEval3DCube(int q, int n, Element<Element_t::CubeEl> const &element = Element<Element_t::CubeEl>());

    BEval3DCube(int q, int n, arma::vec const &coeffVec, Element<Element_t::CubeEl> const &element = Element<Element_t::CubeEl>());

    arma::vec &computeEvaluation();
};

class BEval3DTetra : public BEval<Element_t::TetrahedronEl>
{
  public:
    BEval3DTetra(int q, int n, const Element<Element_t::TetrahedronEl> &element = Element<Element_t::TetrahedronEl>());

    BEval3DTetra(int q, int n, arma::vec const &coeffVec, const Element<Element_t::TetrahedronEl> &element = Element<Element_t::TetrahedronEl>());

    arma::vec &computeEvaluation();
};