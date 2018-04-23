#include <armadillo>
#include "Moments.h"

#ifndef MASSM_H
#define MASSM_H

/*****************************************************************************
 * Bernstein Mass Matrix for 1-dimensional elements                          *
 *****************************************************************************/
class BMass1D : public BMoment1D
{
	int q;				// number of quadrature points ( recommended: 2*(n+1) )
	int n;				// polynomial order
	int lenMass;		// length of the mass matrix
	arma::mat Matrix;	// mass matrix
	int lenBinomialMat; // length of the binomial matrix
	arma::Mat<int64_t> BinomialMat;  // Pascal Matrix

	// computes a Pascal Matrix with size lenBinomialMat
	void compute_binomials();

  public:
	BMass1D(int q, int n);

	~BMass1D();

	// return the mass matrix values at indexes i and j
	double getMatrixValue(int i, int j) { return Matrix(i, j); }

	// compute the mass matrix
	void compute_matrix();

	// sets the function definition or value before computing the mass matrix
	void compute_matrix(double (*f)(double));
	void compute_matrix(arma::vec Fval);
};

/*****************************************************************************
 * Bernstein Mass Matrix for triangular elements (2-dimensional)             *
 *****************************************************************************/
class BMass2DTri : public BMoment2DTri
{
	int q;				// number of quadrature points ( recommended: 2*(n+1) )
	int n;				// polynomial order
	int lenMass;		// length of the matrix
	arma::mat Matrix;	// mass matrix
	int lenBinomialMat; // length of the binomials matrix
	int **BinomialMat;  // computes a Pascal Matrix with size lenBinomialMat

	//alloc matrix linearly
	arma::mat create_matrix();

	void delete_matrix(arma::mat matrix);

	// alloc binomial matrix
	int **create_binomialMat();

	void delete_binomialMat(int **binomialMat);

	// precompute the binomial cofficients necessary
	void compute_binomials();

  public:
	BMass2DTri(int q, int n);

	BMass2DTri(int q, int n, double T[][2]);

	~BMass2DTri();

	double getMatrixValue(int i, int j) { return Matrix(i, j); }

	// compute the mass matrix
	void compute_matrix();

	// sets the function definition or value before computing the mass matrix
	void compute_matrix(double (*f)(double, double));
	void compute_matrix(double *Fval);
};

/*****************************************************************************
 * Bernstein Mass Matrix for quadrilateral elements (2-dimensional)          *
 *****************************************************************************/
class BMass2DQuad : public BMoment2DQuad
{
	int q;				// number of quadrature points ( recommended: 2*(n+1) )
	int n;				// polynomial order
	int lenMass;		// length of the matrix
	arma::mat Matrix;	// mass matrix
	int lenBinomialMat; // length of the binomials matrix
	int **BinomialMat;  // Pascal Matrix, BinomialMat[i, j] == binomial(i+j, i);

	//alloc matrix linearly
	arma::mat create_matrix();

	void delete_matrix(arma::mat matrix);

	// alloc binomial matrix
	int **create_binomialMat();

	void delete_binomialMat(int **binomialMat);

	// computes a Pascal Matrix with size lenBinomialMat
	void compute_binomials();

  public:
	BMass2DQuad(int q, int n);

	BMass2DQuad(int q, int n, double T[][2]);

	~BMass2DQuad();

	double getMatrixValue(int i, int j) { return Matrix(i, j); }

	void compute_matrix();

	void compute_matrix(double (*f)(double, double));

	void compute_matrix(double *Fval);
};

class BMass3D : public BMoment3D
{
};

#endif