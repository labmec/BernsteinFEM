#ifndef ELEM_H
#define ELEM_H

#include "MassM.h"
#include "StiffM.h"

class BElement1D
{
  int q;                // number of element quadrature points
  int n;                // element polynomial basis order
  int length;           // number of element points
  double *BBVector;     // Bernstein-Bézier polynomial basis coefficients
  double *QuadVector;   // BB polynomial evaluated at quadrature points
  double *MassFval;     // mass matrix function values -- stil think it might be unecessary
  //double *ConvecFval;   // convective matrix function values
  double *StiffFval;    // stiffness matrix values
  BMass1D MassMat;     // mass matrix object
  //BConvec1D *ConvecMat; // convective matrix object
  BStiff1D StiffMat;   // stiffness matrix object
  double **ElMat;       // element matrix = mass + convec + stiff
  BMoment1D LoadVec;   // load vector in moments object

  double **create_el_mat();

  void delete_el_mat(double **);

public:
  BElement1D();
  BElement1D(int q, int n);
  ~BElement1D();

  static int sysLen(int n, int nElem) { return n * nElem + 1; }

  int len() { return length; }

  //double** getElementMatrix() { return ElMat; }

  double getMatrixValue (int i, int j) { if (i < length && j < length) return ElMat[i][j]; else return 0.0;}

  double getLoadVector (int i) { if (i < length) return LoadVec.get_bmoment(i); else return 0.0; }

  double getQuadVector(int i) { if (i < q) return QuadVector[i]; else return 0.0; }

  // sets element interval
  void setInterval (double a, double b);

  // sets the BBVector
  void setBBVector (double* vec) { for (int i = 0; i < length; i++) BBVector[i] = vec[i]; }

  // sets functions by function value
  void setMassFunction (double* Fval) { MassMat.setFunction(Fval); }
  //void setConvecFunction (double* Fval);
  void setStiffFunction (double* Fval) {StiffMat.setFunction(Fval); }
  void setLoadFunction (double* Fval) { LoadVec.setFunction(Fval); }
  
  // evaluate the BBVector at the quadrature points
  double* evaluate();

  // computes all matrices and stores the sum into ElMat
  void makeSystem();
};

class BElement2DTri
{
  int q;                   // number of element quadrature points
  int n;                   // element polynomial basis order
  int length;              // number of element points
  double *BBVector;        // Bernstein-Bézier polynomial basis coefficients
  double *QuadVector;      // BB polynomial evaluated at quadrature points
  double *MassFval;        // mass matrix function values
  //double *ConvecFval;      // convective matrix function values
  double *StiffFval;       // stiffness matrix values
  BMass2DTri MassMat;     // mass matrix object
  //BConvec2DTri *ConvecMat; // convective matrix object
  BStiff2DTri StiffMat;   // stiffness matrix object
  double **ElMat;          // element matrix = mass + convec + stiff
  BMoment2DTri LoadVec;   // load vector in moments object

  double **create_el_mat();

  void delete_el_mat(double **);

public:
  BElement2DTri();
  BElement2DTri(int q, int n);
  ~BElement2DTri();

  //double** getElementMatrix() { return ElMat; }

  double getMatrixValue (int i, int j) { if (i < length && j < length) return ElMat[i][j]; else return 0.0;}

  double getLoadVector(int i) { if (i < length) return LoadVec.get_bmoment(i); else return 0.0; }

  double getLoadVector(int a1, int a2) { return LoadVec.get_bmoment(a1, a2); }

  double getQuadVector(int i) { if (i < q) return QuadVector[i]; else return 0.0;}

  // sets element triangle
  void setTriangle (double v1[2], double v2[2], double v3[2]);

  // sets the BBVector
  void setBBVector (double* vec) { for (int i = 0; i < length; i++) BBVector[i] = vec[i]; }

  // sets functions by function value
  void setMassFunction (double* Fval) { MassMat.setFunction(Fval); }
  //void setConvecFunction (double* Fval);
  void setStiffFunction (double* Fval) {StiffMat.setFunction(Fval); }
  void setLoadFunction (double* Fval) { LoadVec.setFunction(Fval); }

  // evaluate the BBVector at the quadrature points
  double* evaluate();

  // computes all matrices and stores the sum into ElMat
  void makeSystem();
};

class BElement2DQuad
{
  int q;                    // number of element quadrature points
  int n;                    // element polynomial basis order
  int length;               // number of element points
  double *BBVector;         // Bernstein-Bézier polynomial basis coefficients
  double *QuadVector;       // BB polynomial evaluated at quadrature points
  double *MassFval;         // mass matrix function values
  //double *ConvecFval;       // convective matrix function values
  double *StiffFval;        // stiffness matrix values
  BMass2DQuad MassMat;     // mass matrix object
  //BConvec2DQuad *ConvecMat; // convective matrix object
  BStiff2DQuad StiffMat;   // stiffness matrix object
  double **ElMat;           // element matrix = mass + convec + stiff
  BMoment2DQuad LoadVec;   // load vector in moments object

  double **create_el_mat();

  void delete_el_mat(double **);

public:
  BElement2DQuad();
  BElement2DQuad(int q, int n);
  ~BElement2DQuad();

  //double** getElementMatrix() { return ElMat; }

  double getMatrixValue (int i, int j) { if (i < length && j < length) return ElMat[i][j]; else return 0.0;}

  double getLoadVector(int i) { if (i < length) return LoadVec.get_bmoment(i); else return 0.0; }

  double getLoadVector(int a1, int a2) { return LoadVec.get_bmoment(a1, a2); }

  double getQuadVector(int i) { if (i < q) return QuadVector[i]; else return 0.0;}

  // sets element quadrilateral
  void setQuad (double v1[2], double v2[2], double v3[2], double v4[2]);

  // sets the BBVector
  void setBBVector (double* vec) { for (int i = 0; i < length; i++) BBVector[i] = vec[i]; }

  // sets functions by function value
  void setMassFunction (double* Fval) { MassMat.setFunction(Fval); }
  //void setConvecFunction (double* Fval);
  void setStiffFunction (double* Fval) {StiffMat.setFunction(Fval); }
  void setLoadFunction (double* Fval) { LoadVec.setFunction(Fval); }
  
  // evaluate the BBVector at the quadrature points
  double* evaluate();

  // computes all matrices and stores the sum into ElMat
  void makeSystem();
};

#endif