#include "MassM.h"

class BMass1D : public BMoment1D
{
	int Mq;
	int Mn;
	double **Matrix;

	void create_matrix()
    {
        double *aux = new double[(n+1) * (n+1)];
        Matrix = new double* [n+1];
        for (int i = 0; i <= n; aux += (n+1), i++)
            Matrix[i] = aux;
    }
public:
	BMass1D (int q, int n) : BMoment1D (q, 2 * n)
	{

    }
	
	~BMass1D()
	{
        delete Matrix[0];
        delete Matrix;
    }


};