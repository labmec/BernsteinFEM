#ifndef DERIVATIVES_H
#define DERIVATIVES_H

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
  QuadDerivative(int q, int n);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  arma::mat getIntegrationPoints();

  void setFunction(const arma::vec &Fval);

  void compute_matrix(const arma::vec &Fval)
  {
    setFunction(Fval);
    compute_matrix();
  }

  virtual void compute_matrix() = 0;
};

class dXi_dXi : public QuadDerivative // does not inherit moment
{
public:
  dXi_dXi(int q, int n);

  // compute matrix coefficients
  void compute_matrix();
};

class dEta_dEta : public QuadDerivative // does not inherit moment
{
public:
  dEta_dEta(int q, int n);

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
  arma::mat Fval;                 // function values at quadrature nodes
  arma::mat vertices;             // Triangle vertices

public:
  TriangleDerivative(int q, int n);

  void setTriangle(const arma::mat &vertices);

  void setTriangle(double v1[], double v2[], double v3[]);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  void setFunction(const arma::mat &Fval);

  void compute_matrix(const arma::mat &Fval)
  {
    setFunction(Fval);
    compute_matrix();
  }

  // compute matrix coefficients
  virtual void compute_matrix() = 0;
};

class dXi_dXi : public TriangleDerivative
{

public:
  dXi_dXi(int q, int n);

  void compute_matrix();
};

class dEta_dEta : TriangleDerivative
{
public:
  dEta_dEta(int q, int n);

  void compute_matrix();
};

class dXi_dEta : TriangleDerivative
{
public:
  dXi_dEta(int q, int n);

  void compute_matrix();
};
} // namespace TriD

#endif