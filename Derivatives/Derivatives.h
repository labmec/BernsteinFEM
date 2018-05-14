#include <armadillo>
#include <stdint.h>
#include "Moments.h"

namespace QuadD
{
class dXi_dXi // does not inherit moment
{
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::mat Fval;                 // function values at quadrature nodes

public:
  dXi_dXi(int q, int n);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  void setFunction(arma::mat Fval);

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval);
};

class dEta_dEta // dows not inherit moment
{
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::mat Fval;                 // function values at quadrature nodes

public:
  dEta_dEta(int q, int n);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  void setFunction(arma::mat Fval);

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval);
};

class dXi_dEta : public BMoment2DQuad
{
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coeffients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix

public:
  dXi_dEta(int q, int n);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval)
  {
    setFunction(Fval);
    compute_matrix();
  }
};
} // namespace QuadD

namespace TriD
{
class dXi_dXi
{
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::mat Fval;                 // function values at quadrature nodes

public:
  dXi_dXi(int q, int n);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  void setFunction(arma::mat Fval);

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval);
};

class dEta_dEta
{
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::mat Fval;                 // function values at quadrature nodes

public:
  dEta_dEta(int q, int n);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  void setFunction(arma::mat Fval);

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval);
};

class dXi_dEta
{
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coefficients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix
  arma::mat Fval;                 // function values at quadrature nodes

public:
  dXi_dEta(int q, int n);

  int Len() { return len; }

  arma::mat getMatrix() { return Matrix; }

  double getMatrixValue(int i, int j) { return Matrix(i, j); }

  void setFunction(arma::mat Fval);

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval);
};
} // namespace TriD