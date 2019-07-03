#pragma once

#include "Elements.h"

// abstract class
template <Element_t EL>
class BEval
{
  protected:
    uint q;                 // number of integration points
    uint n;                 // polynomial order
    arma::vec BBVec;        // Vector holding coefficients
    uint bbvec_len;         // length BBVec should have
    arma::vec eval;         // vector holding the polynomial evaluated in the integration points of the element
    Element<EL> element;    // element of computation
    bool evaluated = false; // indicates whether evaluation was already evaluated or not

  public:
    BEval(uint q, uint n, Element<EL> const &element = Element<EL>());

    BEval(uint q, uint n, arma::vec const &coeffVec, Element<EL> const &element = Element<EL>());

    // getters

    uint getNumIntegrationPoints();

    uint getPOrder();

    arma::vec &getCoefficientsVec();

    arma::vec &getEval();

    Element<EL> &getElement();

    // setters

    void setCoefficientsVec(arma::vec const &vec);

    // sets this object's with values of 0
    void Zero();

    // evaluates the resulting polynomial at the integration points
    // and returns the resulting vector
    virtual arma::vec &computeEvaluation() = 0;

    virtual double L2Norm() = 0;
};

class BEval1D : public BEval<Element_t::LinearEl>
{
  public:
    BEval1D(uint q, uint n, Element<Element_t::LinearEl> const &element = Element<Element_t::LinearEl>());

    BEval1D(uint q, uint n, arma::vec const &coeffVec, Element<Element_t::LinearEl> const &element = Element<Element_t::LinearEl>());

    arma::vec &computeEvaluation();

    double L2Norm();
};

class BEval2DTri : public BEval<Element_t::TriangularEl>
{
    arma::vec eval_inter;

  public:
    BEval2DTri(uint q, uint n, Element<Element_t::TriangularEl> const &element = Element<Element_t::TriangularEl>());

    BEval2DTri(uint q, uint n, arma::vec const &coeffVec, Element<Element_t::TriangularEl> const &element = Element<Element_t::TriangularEl>());

    arma::vec &computeEvaluation();

    double L2Norm();
};

class BEval2DQuad : public BEval<Element_t::QuadrilateralEl>
{
    arma::vec eval_inter;

  public:
    BEval2DQuad(uint q, uint n, Element<Element_t::QuadrilateralEl> const &element = Element<Element_t::QuadrilateralEl>());

    BEval2DQuad(uint q, uint n, arma::vec const &coeffVec, Element<Element_t::QuadrilateralEl> const &element = Element<Element_t::QuadrilateralEl>());

    arma::vec &computeEvaluation();

    double L2Norm();
};

class BEval3DCube : public BEval<Element_t::CubeEl>
{
  public:
    BEval3DCube(uint q, uint n, Element<Element_t::CubeEl> const &element = Element<Element_t::CubeEl>());

    BEval3DCube(uint q, uint n, arma::vec const &coeffVec, Element<Element_t::CubeEl> const &element = Element<Element_t::CubeEl>());

    arma::vec &computeEvaluation();

    double L2Norm();
};

class BEval3DTetra : public BEval<Element_t::TetrahedronEl>
{
  public:
    BEval3DTetra(uint q, uint n, const Element<Element_t::TetrahedronEl> &element = Element<Element_t::TetrahedronEl>());

    BEval3DTetra(uint q, uint n, arma::vec const &coeffVec, const Element<Element_t::TetrahedronEl> &element = Element<Element_t::TetrahedronEl>());

    arma::vec &computeEvaluation();

    double L2Norm();
};