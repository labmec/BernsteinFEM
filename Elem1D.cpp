#include <cmath>

#include "Elem.h"
#include "JacobiGaussNodes.h"

BElement1D::BElement1D()
{
    // get n and q from user, then do the same as the one below
}

BElement1D::BElement1D(int q, int n)
{
    this->q = q;
    this->n = n;

    length = n + 1;

    MassMat = new BMass1D(q, n);
    StiffMat = new BStiff1D(q, n);
    LoadVec = new BMoment1D(q, n);

    QuadVector = new double[q];
    BBVector = new double[length];
    MassFval = new double[length];
    StiffFval = new double[length];

    ElMat = create_el_mat();
}

BElement1D::~BElement1D()
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

double *BElement1D::evaluate()
{
    // initialize with 0's
    for (int i = 0; i < length; i++)
        QuadVector[i] = 0.0;
    
    // evaluate routine, utilizes the recurrence relation of Bernstein Polynomials
    for (int i = 0; i < length; i++)
    {
        double xi = legendre[1][q - 2][i]; // see JacobiGaussNodes.h for reference
        double s = 1 - xi;
        double r = xi / s; // recurrence relation constant
        double w = pow(s, n);
        for (int a = 0; a < length; a++)
        {
            QuadVector[i] += w * BBVector[i];
            w *= r * (n / (1. + a)); // treats the recurrence
        }
    }
    return QuadVector;
}

void BElement1D::makeSystem()
{
    MassMat->compute_matrix();
    //ConvecMat->compute_matrix();
    StiffMat->compute_matrix();
    LoadVec->compute_moments();

    //ElMat = MassMat + StiffMat; //still gonna define the operator overloading
}