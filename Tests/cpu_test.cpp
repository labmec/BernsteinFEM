/** cpu_test.cpp
 *
 * This file consists of a CPU time test for the
 * finite element Stiffness Matrix computation
 * in quadrilateral elements using Bernstein polynomials
 */

#include <iostream>
#include <chrono>

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
    ofstream file("results/cpu_test.txt");
    file << "CPU time results of each iteration of calculating the stiffness matrix" << endl <<
        "Calculating using p+1 integration points and constant function" << endl << endl;
    for(int n = 2; n < 80; n++) {
        file << "P = " << n << ":" << endl;

        int q = n + 1;
        auto init_begin = chrono::high_resolution_clock::now();
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

        auto init_end = chrono::high_resolution_clock::now();

        // computes the stiffness matrix 10 times to make the mean time
        for (int i = 0; i < 10; i++) {
            // computes stiffness matrix
            stif.compute_matrix();

        }
        auto compute_end = chrono::high_resolution_clock::now();
        
        if (n <= 40) {
            auto init_time = (init_end - init_begin);
            auto compute_time = (compute_end - init_end);
            file << "\tInitialization time: " << init_time.count() / 1000000. << " ms" << endl;
            file << "\tComputation time: " << compute_time.count() / 10000000. << " ms" << endl;
        } else {
            auto init_time = (init_end - init_begin);
            auto compute_time = (compute_end - init_end);
            file << "\tInitialization time: " << init_time.count() / 1000000. / 1000 << " s" << endl;
            file << "\tComputation time: " << compute_time.count() / 10000000. / 1000 << " s" << endl;
        }

        
        file << endl;
    }
    file.close();
	return 0;
}