#include "Moments.h"

#ifndef MASSM_H
#define MASSM_H

#ifndef BINOM_PRECOMP
// Usual factorial coefficient
int factorial(int n)
{
	if (n > 1)
		return n * factorial(n - 1);
	else
		return 1;
}

// Overloaded factorial coefficient = factorial(n) / factorial(b),
// used for computing binomial coefficient
int factorial(int a, int b)
{
	if (a > b)
		return a * factorial(a - 1, b);
	else
		return 1;
}

// Usual binomial coefficient
int binomial(int a, int b)
{
	return factorial(a, b) / factorial(a - b);
}
#endif

class BMass1D : public BMoment1D
{
	int q;			 // number of quadrature points
	int n;			 // polynomial order
	int lenMass;	 // length of the mass matrix
	double **Matrix; // mass matrix

	// alloc matrix linearly
	double **create_matrix();

	void delete_matrix(double **matrix);

	// precompute the binomial cofficients necessary
	// not yet implemented in this way
	void compute_binomials();

  public:
	BMass1D(int q, int n);

	~BMass1D();

	// return the mass matrix values at indexes i and j
	double getMatrixValue(int i, int j) { return Matrix[i][j]; }

	// compute the mass matrix
	void compute_matrix();

	// sets the function definition or value before computing the mass matrix
	void compute_matrix(double (*f)(double));
	void compute_matrix(double *Fval);
};

class BMass2DTri : public BMoment2DTri
{
	int q;				// number of quadrature points
	int n;				// polynomial order
	int lenMass;		// length of the matrix
	double **Matrix;	// mass matrix
	int lenBinomialMat; // length of the binomials matrix
	int **BinomialMat;  // computes a Pascal Matrix with size lenBinomialMat

	//alloc matrix linearly
	double **create_matrix();

	void delete_matrix(double **matrix);

	// alloc binomial matrix
	int **create_binomialMat();

	void delete_binomialMat(int **binomialMat);

	// precompute the binomial cofficients necessary
	void compute_binomials();

  public:
	BMass2DTri(int q, int n);

	BMass2DTri(int q, int n, double T[][2]);

	~BMass2DTri();

	static int position(int i1, int j1, int i2, int j2, int n)
	{
		return (i1 + i2) * (n + 1) + (j1 + j2);
	}

	// return the mass matrix values at indexes (i1, j1) and (i2, j2) for 2 triangle coordinates
	double getMatrixValue(int i1, int j1, int i2, int j2);

	// compute the mass matrix
	void compute_matrix();

	// sets the function definition or value before computing the mass matrix
	void compute_matrix(double (*f)(double, double));
	void compute_matrix(double *Fval);
};

class BMass2DQuad : public BMoment2DQuad
{
	int q;
	int n;
	double **Matrix;

  public:
	BMass2DQuad(int q, int n);

	BMass2DQuad(int q, int n, double T[][2]);

	~BMass2DQuad();

	void compute_matrix();

	void compute_matrix(double (*f)(double));

	void compute_matrix(double *Fval);
};

class BMass3D : public BMoment3D
{
};

#endif