#include <cmath>

#include "Elem.h"
#include "JacobiGaussNodes.h"

BElement2DQuad::BElement2DQuad()
{
    // get n and q from user, then do the same as the one below
}

BElement2DQuad::BElement2DQuad(int q, int n)
{
    this->q = q;
    this->n = n;

    length = (n + 1) * (n + 1);

    MassMat = new BMass2DQuad(q, n);
    StiffMat = new BStiff2DQuad(q, n);
    LoadVec = new BMoment2DQuad(q, n);

    QuadVector = new double[q * q];
    BBVector = new double[length];
    MassFval = new double[length];
    StiffFval = new double[length];

    ElMat = create_el_mat();
}

BElement2DQuad::~BElement2DQuad()
{
    delete MassMat;
    delete StiffMat;
    delete BBVector;
    delete QuadVector;
    delete MassFval;
    //delete ConvecFval;
    delete StiffFval;
    delete_el_mat(ElMat);
}

double *BElement2DQuad::evaluate()
{
    // initialize QuadVector with 0's
    for (int i = 0; i < q * q; i++)
        QuadVector[i] = 0.0;
    
    // EvalStep routine, same as the 1D with tensor product
    for (int i = 0; i < q; i++)
    {
        double xi1 = legendre[1][q - 2][i]; // see JacobiGaussNodes.h for reference
        double s1 = 1 - xi1;
        double r1 = xi1 / s1;
        double w1 = pow(s1, n);
        for (int j = 0; j < q; j++)
        {
            double xi2 = legendre[1][q - 2][j];
            double s2 = 1 - xi2;
            double r2 = xi2 / s2;
            double w2 = pow(s2, n);

            for (int a = 0; a < n+1; a++)
            {
                for (int b = 0; b < n+1; b++)
                {
                    QuadVector[n * i + j] += w2 * BBVector[n * i + j];
                    w2 *= r2 * (n / (1. + b));
                }
                w1 *= r1 * (n / 1. + a);
            }
        }
    }
    return QuadVector;
}

void BElement2DQuad::makeSystem()
{
    MassMat->compute_matrix();
    //ConvecMat->compute_matrix();
    StiffMat->compute_matrix();
    LoadVec->compute_moments();

    //ElMat = MassMat + StiffMat; //still gonna define the operator overloading
}