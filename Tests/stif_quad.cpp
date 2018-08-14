#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <armadillo>

#include "Derivatives.h"
#include "JacobiGaussNodes.h"

using namespace std;
using namespace QuadD;

double f(double x, double y)
{
    return  x * y + (1 - x * x);
}

int main()
{
    int n = 2;
    int q = 4 * n;

    StiffnessMatrix stif(q, n);

    ofstream file;
    
    file.open("results/stif_quad.txt");

    // sets quadrilateral
    {
        // arma::mat vertices = {{0,0},{1,0},{1,1},{0,1}};
        // arma::mat vertices = {{1, 1}, {5, -2}, {3, 2}, {0, 2}};
        arma::mat vertices = {{1., 1.}, {3., 2.}, {3., 3.}, {1., 3.}};
        stif.setQuadrilateral(vertices);
    }
    
    
    // arma::mat quadPoints = stif.getIntegrationPoints();

    // sets function values
    {
        arma::vec Fval(q * q, arma::fill::ones);
        
        // arma::vec Fval(quadPoints.n_rows, arma::fill::none);
        
        // for (int i = 0; i < Fval.n_rows; i++)
        // {
        //     Fval(i) = f(quadPoints(i, 0), quadPoints(i, 1));
        // }
        
        stif.setFunction(Fval);
    }

    // computes stiffness matrix
    stif.compute_matrix();

    // prints stiffness matrix to file
    file << "{";
    for (int i = 0; i < (n+1) * (n+1); i++)
    {
        file << "{";
        for (int j = 0; j < (n+1) * (n+1); j++)
            if (j < (n+1)*(n+1) - 1)
                file << stif.getMatrixValue(i, j) << ", ";
            else   
                file << stif.getMatrixValue(i, j);
        if (i < (n+1)*(n+1) - 1)
            file << "},\n";
        else
            file << "}";
    }
    file << "}" << endl;

	return 0;
}