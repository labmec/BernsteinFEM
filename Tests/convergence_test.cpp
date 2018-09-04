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
    int n = 8;

    ofstream file("results/convergence_test.txt");

    file << "{";
    for(int q = 2; q < 40; q += 1) {
        StiffnessMatrix stif(q, n);

        // sets quadrilateral
        {
            arma::mat vertices = {{0,0},{1,0},{1,1},{0,1}};
            // arma::mat vertices = {{1, 1}, {5, -2}, {3, 2}, {0, 2}};
            // arma::mat vertices = {{1., 1.}, {3., 2.}, {3., 3.}, {1., 3.}};
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

        stif.print_mathematica(file);
        if (q < 39) {
            file << ",";
        }
        file << endl;
    }
    file << "}";
    file.close();
	return 0;
}