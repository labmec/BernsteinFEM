#pragma once

#include <armadillo>
#include <stdint.h>
#include "Moments.h"

namespace QuadD
{
// abstract class -- do not instantiate
class QuadDerivative
{
  void compute_binomials(arma::Mat<int64_t> &BinomialMat, int lenBinom);

protected:
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::vec Fval;                 // function values at quadrature nodes

public:
  // constructor
  QuadDerivative(int q, int n);

  // returns the length of the square Matrix
  int Len() { return len; }

  // returns a reference to the Matrix
  const arma::mat &getMatrix() { return Matrix; }

  // returns the Matrix[i][j] value
  double getMatrixValue(int i, int j) { return Matrix(i, j); }

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
  dXi_dXi(int q, int n);

  void setFunction(const arma::vec &Fval)
  {
    BMoment2DQuad::setFunction(Fval);
  }

  // compute matrix coefficients
  void compute_matrix();
};

class dEta_dEta : private BMoment2DQuad, public QuadDerivative // does not inherit moment
{
public:
  dEta_dEta(int q, int n);

  void setFunction(const arma::vec &Fval)
  {
    BMoment2DQuad::setFunction(Fval);
  }

  // compute matrix coefficients
  void compute_matrix();
};

class dXi_dEta : private BMoment2DQuad, public QuadDerivative
{
public:
  dXi_dEta(int q, int n);

  void setFunction(const arma::vec &Fval)
  {
    BMoment2DQuad::setFunction(Fval);
  }

  // compute matrix coefficients
  void compute_matrix();
};

class StiffnessMatrix : public QuadDerivative
{
  dXi_dXi Xi_Xi;
  dXi_dEta Xi_Eta;
  dEta_dEta Eta_Eta;

  arma::mat vertices;

public:
  StiffnessMatrix(int q, int n);

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

namespace TriD
{

// abstract class -- do not instantiate
class TriangleDerivative
{
  void compute_binomials(arma::Mat<int64_t> &BinomialMat, int lenBinom);

protected:
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::vec Fval;                 // function values at quadrature nodes
  arma::mat vertices;             // Triangle vertices

public:
  TriangleDerivative(int q, int n);

  void setTriangle(const arma::mat &vertices);

  void setTriangle(double v1[], double v2[], double v3[]);

  double area();

  int Len() { return len; }

  void Zero() { Matrix.zeros(); }

  const arma::mat &getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  void setFunction(const arma::mat &Fval);

  void compute_matrix(const arma::mat &Fval)
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

  // compute matrix coefficients
  virtual void compute_matrix() = 0;

  const double &operator() (int i, int j);
};

class dXi_dXi : public TriangleDerivative
{

public:
  dXi_dXi(int q, int n);

  void compute_matrix();
};

class dEta_dEta : public TriangleDerivative
{
public:
  dEta_dEta(int q, int n);

  void compute_matrix();
};

class dXi_dEta : public TriangleDerivative
{
public:
  dXi_dEta(int q, int n);

  void compute_matrix();
};

class StiffnessMatrix : public TriangleDerivative
{
  dXi_dXi Xi_Xi;
  dXi_dEta Xi_Eta;
  dEta_dEta Eta_Eta;

  arma::mat vertices;

  // makes Duffy Transform with respect to the Triangle in this object
  // using points X, returns the values in points Xi
  void duffyTransform (const arma::vec &X, arma::vec &Xi);
  
  // makes Duffy Transform with respect to the Triangle in this object
  // using parameter b1,b2,b3, returns the values in points Xi
  void duffyTransform (double b1, double b2, double b3, arma::vec &Xi);

public:
  StiffnessMatrix(int q, int n);

  void Zero()
  {
    Xi_Xi.Zero();
    Xi_Eta.Zero();
    Eta_Eta.Zero();
    TriangleDerivative::Zero();
  }

  void setFunction(const arma::vec &Fval);

  void setTriangle(const arma::mat &vertices);

  arma::mat getIntegrationPoints();

  // compute matrix coefficients
  void compute_matrix();
};
} // namespace TriD
