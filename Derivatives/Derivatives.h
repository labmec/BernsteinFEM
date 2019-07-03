#pragma once

#include <stdint.h>
#include "Moments.h"

namespace QuadD
{
// abstract class -- do not instantiate
class QuadDerivative
{
  void compute_binomials(arma::Mat<int64_t> &BinomialMat, uint lenBinom);

protected:
  uint q;                          // number of quadrature nodes
  uint n;                          // polynomial order
  uint len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  uint lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::vec Fval;                 // function values at quadrature nodes

public:
  // constructor
  QuadDerivative(uint q, uint n);

  // returns the length of the square Matrix
  uint Len() { return len; }

  // returns a reference to the Matrix
  const arma::mat &getMatrix() { return Matrix; }

  // returns the Matrix[i][j] value
  double getMatrixValue(uint i, uint j) { return Matrix(i, j); }

  // zeroes the Matrix
  void Zero() { Matrix.zeros(); }

  // returns the integration points to the quadrilateral
  arma::mat getIntegrationPoints();

  void setFunction(const arma::vec &Fval);

  void compute_matrix(const arma::vec &Fval)
  {
    setFunction(Fval);
    compute_matrix();
  }

  // prints the Matrix in format compatible with Mathematica to the specified stream
  // default stream: cout
  void print_mathematica(std::ostream &stream = std::cout);

  // prints the Matrix in readable format to the specified stream
  // default stream: cout
  void print(std::ostream &stream = std::cout);

  virtual void compute_matrix() = 0;
};

class dXi_dXi : private BMoment2DQuad, public QuadDerivative // does not inherit moment
{
public:
  dXi_dXi(uint q, uint n);

  void setFunction(const arma::vec &Fval)
  {
    BMoment2DQuad::setFunctionValues(Fval);
  }

  // compute matrix coefficients
  void compute_matrix();
};

class dEta_dEta : private BMoment2DQuad, public QuadDerivative // does not inherit moment
{
public:
  dEta_dEta(uint q, uint n);

  void setFunction(const arma::vec &Fval)
  {
    BMoment2DQuad::setFunctionValues(Fval);
  }

  // compute matrix coefficients
  void compute_matrix();
};

class dXi_dEta : private BMoment2DQuad, public QuadDerivative
{
public:
  dXi_dEta(uint q, uint n);

  void setFunction(const arma::vec &Fval)
  {
    BMoment2DQuad::setFunctionValues(Fval);
  }

  // compute matrix coefficients
  void compute_matrix();
};

class StiffnessMatrix : public QuadDerivative
{
  dXi_dXi Xi_Xi;
  dXi_dEta Xi_Eta;
  dEta_dEta Eta_Eta;

  Element<Element_t::QuadrilateralEl> element;

public:
  StiffnessMatrix(uint q, uint n, const Element<Element_t::QuadrilateralEl> &element = Element<Element_t::QuadrilateralEl>());

  void Zero()
  {
    Xi_Xi.Zero();
    Xi_Eta.Zero();
    Eta_Eta.Zero();
    QuadDerivative::Zero();
  }

  void setFunction(const arma::vec &Fval);

  void setQuadrilateral(const arma::mat &vertices);

  arma::mat getIntegrationPoints();

  // compute matrix coefficients
  void compute_matrix();
};
} // namespace QuadD