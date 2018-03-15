#include "Moments.h"

#ifndef MASSM_H
#define MASSM_H

/*to implement yet
  About computing the Element-Mass matrix
*/
class BMass1D : public BMoment1D
{
	int Mq;
	int Mn;
	double **Matrix;

	void create_matrix();
public:
	BMass1D() : BMoment1D()
	{	}

	BMass1D (int q, int n) : BMoment1D (q, 2 * n)
	{	}
	
	~BMass1D();
};

class BMass2DTri : public BMoment2DTri
{
public:
    // check about the 'q' value
    BMass2DTri(int q, int n) : BMoment2DTri(q, 2 * n)
    {	}

    BMass2DTri(int q, int n, double T[][2]) : BMoment2DTri(q, 2 * n, T)
	{	}
};

class BMass2DQuad : public BMoment2DQuad
{

};

class BMass3D : public BMoment3D
{

};

#endif