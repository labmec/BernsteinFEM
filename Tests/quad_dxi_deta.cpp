/** quad_dxi_deta.cpp
 * 
 * This file makes a test for the computation
 * of the dXi_dEta object matrix in quadrilateral
 * elements
 */

#include "Derivatives.h"
#include <armadillo>
#include <iostream>

using namespace QuadD;
using namespace arma;
using namespace std;

int main() {
    int n = 2;
    int q = 2 * n;
    dXi_dEta X(q, n);
    int len = X.Len();
    ofstream file;

    file.open("results/dxi_deta.txt");
    
    mat funct_val(q * q, 1, fill::ones);

    X.setFunction(funct_val);

    X.compute_matrix();

    file << "{";
    for (int i = 0; i < len; i++)
    {
        file << "{";
        for (int j = 0; j < len - 1; j++)
        {
            file << X.getMatrixValue(i, j) << ",\t";
        }
        if (i < len - 1) {
            file << X.getMatrixValue(i, len - 1) << "},\n";
        } else {
            file << X.getMatrixValue(i, len - 1) << "}";
        }
    }
    file << "}";

    return 0;
}