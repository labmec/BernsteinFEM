#include <iostream>

#define USE_MATH_DEFINES
#include <cmath>

#include <armadillo>

#include "Elem1D.cpp"
#include "JacobiGaussNodes.h"

using namespace std;

int main()
{
    int n = 4;
    int q = 2 * n;

    BElement1D elem(q, n);

    // set element interval
    {
        double a = 0.0;
        double b = 1.0;
        elem.setInterval(a, b);
    }

    // set functions values (mass, stiffness and load vec) all as 1
    {
        arma::vec Fval(q, arma::fill::ones);
        
        elem.setLoadFunction(Fval);
        elem.setMassFunction(Fval);
        elem.setStiffFunction(Fval);
    }

    // compute element mass matrix, stiffness matrix and load vector
    elem.makeSystem();

    // print element matrix == MassMat + StiffMAt
    cout << "Element's matrix:" << endl << endl;
    for (int i = 0; i < elem.len(); i++)
    {
        for (int j = 0; j < elem.len(); j++)
            cout << scientific << elem.getMatrixValue(i, j) << " ,  ";
        cout << endl;
    }
    
    // print load vector
    cout << endl << "Element's Load Vector" << endl << endl;
    for (int i = 0; i < elem.len(); i++)
        cout << "LVec[" << i << "]: " << scientific << elem.getLoadVector(i) << endl;
    
    return 0;
}