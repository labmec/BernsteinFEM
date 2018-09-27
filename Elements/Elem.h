#include <armadillo>
#include "MassM.h"
#include "StiffM.h"

#ifndef ELEM_H
#define ELEM_H

class BElement1D
{
  int q;                // number of element quadrature points
  int n;                // element polynomial basis order
  int length;           // number of element points
  arma::vec BBVector;   // Bernstein-Bézier polynomial basis coefficients
  arma::vec QuadVector; // BB polynomial evaluated at quadrature points
  arma::vec MassFval;   // mass matrix function values -- stil think it might be unecessary
  //arma::vec ConvecFval; // convective matrix function values
  arma::vec StiffFval;  // stiffness matrix values
  BMass1D MassMat;      // mass matrix object
  //BConvec1D *ConvecMat; // convective matrix object
  BStiff1D StiffMat;    // stiffness matrix object
  arma::mat ElMat;      // element matrix = mass + convec + stiff
  BMoment1D LoadVec;    // load vector in moments object

public:
  BElement1D();
  BElement1D(int q, int n);
  ~BElement1D();

  static int sysLen(int n, int nElem)
  {
    return n * nElem + 1;
  }

  //double** getElementMatrix() { return ElMat; }

  int len() { return length; }

  double getMatrixValue(int i, int j)
  {
    return ElMat(i, j);
  }

  double getLoadVector(int i)
  {
    if (i < length)
      return LoadVec.get_bmoment(i);
    else
      return 0.0;
  }

  double getQuadVector(int i)
  {
    if (i < q)
      return QuadVector[i];
    else
      return 0.0;
  }

  // sets element interval
  void setInterval(double a, double b);

  // sets the BBVector
  void setBBVector(double *vec)
  {
    for (int i = 0; i < length; i++)
      BBVector(i) = vec[i];
  }

  void setBBVector(arma::vec vec)
  {
    for (int i = 0; i < length; i++)
      BBVector(i) = vec(i);
  }

  // sets functions by function value
  void setMassFunction(arma::vec Fval) { MassMat.setFunction(Fval); }
  //void setConvecFunction (double* Fval);
  void setStiffFunction(arma::vec Fval) { StiffMat.setFunction(Fval); }
  void setLoadFunction(arma::vec Fval) { LoadVec.setFunction(Fval); }

  // evaluate the BBVector at the quadrature points
  arma::vec evaluate();

  // computes all matrices and stores the sum into ElMat
  void makeSystem();
};

class BElement2DTri
{
  int q;                   // number of element quadrature points
  int n;                   // element polynomial basis order
  int length;              // number of element points
  arma::vec BBVector;      // Bernstein-Bézier polynomial basis coefficients
  arma::vec QuadVector;    // BB polynomial evaluated at quadrature points
  arma::vec MassFval;      // mass matrix function values
  //arma::vec ConvecFval;    // convective matrix function values
  arma::vec StiffFval;     // stiffness matrix values
  BMass2DTri MassMat;      // mass matrix object
  //BConvec2DTri *ConvecMat; // convective matrix object
  BStiff2DTri StiffMat;    // stiffness matrix object
  arma::mat ElMat;         // element matrix = mass + convec + stiff
  BMoment2DTri LoadVec;    // load vector in moments object

  arma::mat create_el_mat();

  void delete_el_mat(arma::mat);

public:
  BElement2DTri();
  BElement2DTri(int q, int n);
  ~BElement2DTri();

  //double** getElementMatrix() { return ElMat; }

  int len() { return length; }

  double getMatrixValue(int i, int j)
  {
    return ElMat(i, j);
  }

  double getLoadVector(int i)
  {
    if (i < length)
      return LoadVec.get_bmoment(i);
    else
      return 0.0;
  }

  double getLoadVector(int a1, int a2) { return LoadVec.get_bmoment(a1, a2); }

  double getQuadVector(int i)
  {
    if (i < q)
      return QuadVector[i];
    else
      return 0.0;
  }

  // sets element triangle
  void setTriangle(double v1[2], double v2[2], double v3[2]);

  // sets the BBVector
  void setBBVector(double *vec)
  {
    for (int i = 0; i < length; i++)
      BBVector[i] = vec[i];
  }

  // sets functions by function value
  void setMassFunction(arma::vec Fval) { MassMat.setFunction(Fval); }
  //void setConvecFunction (double* Fval);
  void setStiffFunction(arma::vec Fval) { StiffMat.setFunction(Fval); }
  void setLoadFunction(arma::vec Fval) { LoadVec.setFunction(Fval); }

  // evaluate the BBVector at the quadrature points
  arma::vec evaluate();

  // computes all matrices and stores the sum into ElMat
  void makeSystem();
};

class BElement2DQuad
{
  int q;                    // number of element quadrature points
  int n;                    // element polynomial basis order
  int length;               // number of element points
  arma::vec BBVector;       // Bernstein-Bézier polynomial basis coefficients
  arma::vec QuadVector;     // BB polynomial evaluated at quadrature points
  arma::vec MassFval;       // mass matrix function values
  //arma::vec ConvecFval;     // convective matrix function values
  arma::vec StiffFval;      // stiffness matrix values
  BMass2DQuad MassMat;      // mass matrix object
  //BConvec2DQuad *ConvecMat; // convective matrix object
  BStiff2DQuad StiffMat;    // stiffness matrix object
  arma::mat ElMat;          // element matrix = mass + convec + stiff
  BMoment2DQuad LoadVec;    // load vector in moments object

  arma::mat create_el_mat();

  void delete_el_mat(arma::mat);

public:
  BElement2DQuad();
  BElement2DQuad(int q, int n);
  ~BElement2DQuad();

  //double** getElementMatrix() { return ElMat; }

  double getMatrixValue(int i, int j)
  {
    return ElMat(i, j);
  }

  double getLoadVector(int i)
  {
    if (i < length)
      return LoadVec.get_bmoment(i);
    else
      return 0.0;
  }

  double getLoadVector(int a1, int a2) { return LoadVec.get_bmoment(a1, a2); }

  double getQuadVector(int i)
  {
    if (i < q)
      return QuadVector[i];
    else
      return 0.0;
  }

  // sets element quadrilateral
  void setQuad(double v1[2], double v2[2], double v3[2], double v4[2]);

  // sets the BBVector
  void setBBVector(double *vec)
  {
    for (int i = 0; i < length; i++)
      BBVector(i) = vec[i];
  }

  void setBBVector(arma::vec vec)
  {
    for (int i = 0; i < length; i++)
      BBVector(i) = vec(i);
  }

  // sets functions by function value
  void setMassFunction(arma::vec Fval) { MassMat.setFunction(Fval); }
  //void setConvecFunction (double* Fval);
  void setStiffFunction(arma::vec Fval) { StiffMat.setFunction(Fval); }
  void setLoadFunction(arma::vec Fval) { LoadVec.setFunction(Fval); }

  // evaluate the BBVector at the quadrature points
  arma::vec evaluate();

  // computes all matrices and stores the sum into ElMat
  void makeSystem();
};

#endif