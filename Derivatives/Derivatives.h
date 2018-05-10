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

  void setFunction(arma::mat Fval);

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval);
}

class dEta_dEta // dows not inherit moment
{

}

class dXi_dEta : public BMoment2DQuad
{
  int q;                          // number of quadrature nodes
  int n;                          // polynomial order
  int len;                        // length of the matrix
  arma::mat Matrix;               // matrix containing coeffients
  int lenBinom;                   // length of the binomial (Pascal) matrix
  arma::Mat<int64_t> BinomialMat; // Pascal matrix

  // computes object's Pascal Matrix
  void compute_binomials();

public:
  xXi_dEta(int q, int n);

  // compute matrix coefficients
  void compute_matrix();
  void compute_matrix(arma::mat Fval)
  {
    setFunction(Fval);
    compute_matrix();
  }
}
} // namespace QuadD

namespace TriD
{
class dXi_dXi
{
}

class dEta_dEta
{

}

class dXi_dEta
{
}
} // namespace TriD