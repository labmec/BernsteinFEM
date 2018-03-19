#include "Moments.h"

#ifndef MASSM_H
#define MASSM_H

/*to implement yet
  About computing the Element-Mass matrix
*/
class BMass1D : public BMoment1D
{
	int q;
	int n;
	int lenMass;
	double **Matrix;
	int *BinomialMat;

	// alloc matrix linearly
	double** create_matrix();

	void compute_binomials();

	// Usual factorial coefficient
	int factorial (int n);
	
	// Overloaded factorial coefficient = factorial(n) / factorial(b),
	// used for computing binomial coefficient
	int factorial (int a, int b);

	// Usual binomial coefficient
	double binomial (int a, int b);

public:
	BMass1D (int q, int n);
	
	~BMass1D ();

	// return the mass matrix values at indexes i and j
	double getMatrixValue(int i, int j) { return Matrix[i][j]; }

	// compute the mass matrix
	void compute_matrix();

	// sets the function definition or value before computing the mass matrix
	void compute_matrix (double (*f) (double));
	void compute_matrix (double *Fval);
};

class BMass2DTri : public BMoment2DTri
{
	int q;
	int n;
	double **Matrix;

public:
    BMass2DTri (int q, int n);

    BMass2DTri (int q, int n, double T[][2]);

	~BMass2DTri ();

	void compute_matrix ();

	void compute_matrix (double (*f) (double));

	void compute_matrix (double *Fval);
};

class BMass2DQuad : public BMoment2DQuad
{
	int q;
	int n;
	double **Matrix;

public:
	BMass2DQuad (int q, int n);

    BMass2DQuad (int q, int n, double T[][2]);

	~BMass2DQuad ();

	void compute_matrix ();

	void compute_matrix (double (*f) (double));

	void compute_matrix (double *Fval);
};

class BMass3D : public BMoment3D
{

};

#endif