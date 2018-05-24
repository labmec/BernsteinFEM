#include "Derivatives.h"
#include <armadillo>
#include <iostream>

using namespace QuadD;
using namespace arma;
using namespace std;

int main() {
    int n = 2;
    int q = 2 * n;
    dXi_dXi X(q, n);
    int len = X.Len();
    
    mat funct_val(q, q, fill::ones);

    X.setFunction(funct_val);

    X.compute_matrix();

    for (int i = 0; i < len; i++)
    {
        for (int j = 0; j < len; j++)
        {
            cout << scientific << X.getMatrixValue(i, j) << ";\t";
        }
        cout << endl;
    }
    return 0;
}