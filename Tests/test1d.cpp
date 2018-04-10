#include <cstdio>
#include <iostream>
#include <vector>

#define USE_MATH_DEFINES
#include <cmath>

#include "Elem.h"
#include "JacobiGaussNodes.h"

using namespace std;

// load function
double loadF (double x)
{
    return M_PI * M_PI * sin(M_PI * x);
}

void assemble(BElement1D, double**, double*, int, int);
double* linearSolve(double**, double*, int);
void printEl(BElement1D);

int main()
{
    int n = 4; // polynomials order
    int q = 2 * n; // quadrature points
    int nElem = 10; // number of elements
    int len = BElement1D::sysLen(n, nElem); //length of the global system

    // boundary variables
    double a = 0.0, b = 1.0;
    double h = (b - a) / nElem;

    // instatiate elements objects (they have most of the methods needed)
    std::vector <BElement1D> Elements (nElem, BElement1D(q, n));

    // sets the intervals for each element
    {
        double ah = a;
        double bh = a + h;
        for (int i = 0; i < nElem; i++)
        {
            Elements[i].setInterval(ah, bh);
        }
    }
    
    // compute functions values at quadrature points
    {
        double **lvec = new double*[nElem];
        for (int k = 0; k < nElem; k++)
            lvec[k] = new double[q];
        double *stiff = new double[6 * q]; //must have 6 times the size
        double *mass = new double[q];
        for (int i = 0; i < q; i++)
        {
            double ah = a;
            for (int k = 0; k < nElem; k++)
            {
                lvec[k][i] = loadF( (legendre_xi(q, i) + ah) * (h/2) );
                ah += h;
            }
            mass[i] = 0.0;
            stiff[i] = stiff[q + i] = stiff[ q+q+i ] = stiff[ q+q+q+i ] = stiff[ q+q+q+q+i ] = stiff[ q+q+q+q+q+i ] = 1.0;
        }
        for (int k = 0; k < nElem; k++)
        {
            Elements[k].setStiffFunction(stiff);
            Elements[k].setMassFunction(mass);
            Elements[k].setLoadFunction(lvec[k]);
        }
        
        for (int k = 0; k < nElem; k++)
            delete lvec[k];
        delete lvec;
        delete stiff;
        delete mass;
    }

    // computes elements matrices
    for (int i = 0; i < nElem; i++)
        Elements[i].makeSystem();

    // print elements matrices
    for(int i = 0; i < nElem; i++)
        printEl(Elements[i]);
    
    return 0;
}

void printEl(BElement1D el)
{
    cout << "Element matrix: " << endl << endl;
    for (int i = 0; i < el.len(); i++)
    {
        for (int j = 0; j < el.len(); j++)
            cout << scientific << el.getMatrixValue(i, j) << "  ";
        cout << endl;
    }
    cout << endl;
}