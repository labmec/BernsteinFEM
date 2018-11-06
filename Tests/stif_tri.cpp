/** stif_quad.cpp
 * 
 * This file makes tests for the computation
 * of the finite element Stiffness Matrix in
 * quadrilateral elements using Bernstein polynomials
 * 
 * Remark: this test is not made using the Stiffness directory
 *  code but using the class defined in the Derivatives directory
 *  furthermore this shall change into the Stiffness directory
 */

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <armadillo>

#include "StiffM.h"
#include "JacobiGaussNodes.h"

using namespace std;


double f(double x, double y)
{
    return  x * y + (1 - x * x);
}

int main()
{
    int n = 2;
    int q = 2 * n;

    BStiff2DTri stif(q, n);

    // sets triangle
    {
        arma::mat vertices = {{0.,0.},{2.,0.},{0.,1.}};

        stif.setTriangle(vertices);
    }
    
    
    arma::mat quadPoints = stif.getIntegrationPoints();

    // sets function values
    {
        // arma::vec Fval(quadPoints.n_rows, arma::fill::ones);
        
        arma::vec Fval(quadPoints.n_rows, arma::fill::none);
        
        for (u_int i = 0; i < Fval.n_rows; i++)
        {
            Fval(i) = f(quadPoints(i, 0), quadPoints(i, 1));
        }
        quadPoints.print("quadPoints:");
        stif.setFunction(Fval);
    }

    // computes stiffness matrix
    stif.compute_matrix();

    ofstream file;

    file.open("results/stif_tri.txt");
    file.precision(5);
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

    file.close();

	return 0;
}