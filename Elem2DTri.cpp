#include <cmath>

#include "Elem.h"
#include "JacobiGaussNodes.h"

BElement2DTri::BElement2DTri()
{
    // get n and q from user, then do the same as the one below
}

BElement2DTri::BElement2DTri(int q, int n)
{
    this->q = q;
    this->n = n;

    //length = ((n + 1) * (n + 2)) / 2;
    length = (n + 1) * (n + 1); // accounts for positioning

    MassMat = new BMass2DTri(q, n);
    StiffMat = new BStiff2DTri(q, n);
    LoadVec = new BMoment2DTri(q, n);

    QuadVector = new double[q * q];
    BBVector = new double[length];
    MassFval = new double[length];
    StiffFval = new double[length];

    ElMat = create_el_mat();
}

BElement2DTri::~BElement2DTri()
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

double *BElement2DTri::evaluate()
{
    int Length = ((n + 1) * (n + 2)) / 2;
    double *QuadInter = new double[q * q];
    for (int i = 0; i < q; i++)
        QuadVector[i] = QuadInter[i] = 0.0;

    // convert first index
    for (int i = 0; i < q; i++)
    {
        double xi = legendre[1][q - 2][i]; // see JacobiGaussNodes.h for reference
        double s = 1 - xi;
        double r = xi / s;

        for (int a1 = 0; a1 < Length; a1++)
        {
            double w = pow(s, n - a1);
            for (int a2 = 0; a2 < n - a1; a2++)
            {
                QuadInter[i] += w * BBVector[BMoment2DTri::position(a1, a2, n)];
                w *= r * (n - a1 - a2) / (1. + a2);
            }
        }
    }

    // convert second index
    for (int i = 0; i < q; i++)
    {
        double xi = jacobi[1][q - 2][i]; // see JacobiGaussNodes.h for reference
        double s = 1 - xi;
        double r = xi / s;
        double w = pow(s, n);
        for (int a1 = 0; a1 < Length; a1++)
        {
            for (int i2 = 0; i2 < q; i2++)
            {
                QuadVector[i] += w * QuadInter[i2];
            }
            w *= r * (n -a1 / (1. + a1));
        }
    }

    delete QuadInter;

    return QuadVector;
}

void BElement2DTri::makeSystem()
{
    MassMat->compute_matrix();
    //ConvecMat->compute_matrix();
    StiffMat->compute_matrix();
    LoadVec->compute_moments();

    //ElMat = MassMat + StiffMat; //still gonna define the operator overloading
}