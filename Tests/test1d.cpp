#include <cstdio>
#include <iostream>
#include <vector>
#include <chrono>

#define USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "Elem.h"
#include "JacobiGaussNodes.h"

#define MAX_ITERATIONS 100000

using namespace std;

// load function
double loadF(double x)
{
    return M_PI * M_PI * sin(M_PI * x);
}

void assemble(BElement1D **, mat, double *, int, int);
double *linearSolve(mat, double *, int);
void printEl(BElement1D, int);

int main()
{
    int n = 8;                              // polynomials order
    int q = 2 * n;                          // quadrature points
    int nElem = 100000;                         // number of elements
    int len = BElement1D::sysLen(n, nElem); //length of the global system

    // boundary variables
    double a = 0.0, b = 1.0;
    double h = (b - a) / nElem;

    // instatiate elements objects (they have most of the methods needed)
    //std::vector <BElement1D> Elements (nElem, BElement1D(q, n));
    // somehow the allocator is not allocating the Matrices

    BElement1D **Elements = new BElement1D *[nElem];

    for (int i = 0; i < nElem; i++)
        Elements[i] = new BElement1D(q, n);

    // sets the intervals for each element
    {
        double ah = a;
        double bh = a + h;
        for (int i = 0; i < nElem; i++)
        {
            Elements[i]->setInterval(ah, bh);
        }
    }

    auto start = chrono::high_resolution_clock::now();
    // compute functions values at quadrature points
    {
        matlvec = new double *[nElem];
        for (int k = 0; k < nElem; k++)
            lvec[k] = new double[q];
        double *stiff = new double[6 * q]; //must have 6 times the size
        double *mass = new double[q];
        for (int i = 0; i < q; i++)
        {
            double ah = a;
            for (int k = 0; k < nElem; k++)
            {
                lvec[k][i] = loadF((legendre_xi(q, i) + ah) * (h / 2));
                ah += h;
            }
            mass[i] = 0.0;
            stiff[i] = stiff[q + i] = stiff[q + q + i] = stiff[q + q + q + i] = stiff[q + q + q + q + i] = stiff[q + q + q + q + q + i] = 1.0;
        }
        for (int k = 0; k < nElem; k++)
        {
            Elements[k]->setStiffFunction(stiff);
            Elements[k]->setMassFunction(mass);
            Elements[k]->setLoadFunction(lvec[k]);
        }

        for (int k = 0; k < nElem; k++)
            delete lvec[k];
        delete lvec;
        delete stiff;
        delete mass;
    }
    auto end = chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = end - start;
    std::cout << "Time to setup all elements: " << diff.count() << " s\n";

    start = chrono::high_resolution_clock::now();
    // computes elements matrices
    for (int i = 0; i < nElem; i++)
        Elements[i]->makeSystem();
    end = chrono::high_resolution_clock::now();

    // print elements matrices
    //for(int i = 0; i < nElem; i++)
    //    printEl(*Elements[i], i);

    diff = end - start;
    std::cout << "Time to compute all elements: " << diff.count() << " s\n";

    matGlobalMatrix = new double *[len];
    for (int i = 0; i < len; i++)
        GlobalMatrix[i] = new double[len];
    double *GlobalVec = new double[len];

    // assemble global linear system
    assemble(Elements, GlobalMatrix, GlobalVec, n, nElem);
/*
    // print global matrix and load vector
    cout << endl << "Global Matrix: " << endl << endl << "{";
    for (int i = 0; i < len; i++)
    {
        cout << "{";
        for (int j = 0; j < len; j++)
            if (j < len - 1)
                cout << GlobalMatrix[i][j] << ", ";
            else
                cout << GlobalMatrix[i][j];
        if (i < len - 1)
            cout << "},";
        else
            cout << "}";
    }
    cout << "}";
    cout << endl << "Load Vec: " << endl << endl << "{";
    for (int i = 0; i < len; i++)
        if (i < len - 1)
            cout << GlobalVec[i] << ", ";
        else
            cout << GlobalVec[i];
    cout << "}" << endl << endl;

    // solve global linear system
    double *GlbBBVector = linearSolve(GlobalMatrix, GlobalVec, len);

    cout << "Resolution, global BBVector: " << endl;
    for (int i = 0; i < len; i++)
        cout << GlbBBVector[i] << " ;  ";
    cout << endl
         << endl;
*/
    return 0;
}

void printEl(BElement1D el, int i)
{
    cout << "Element matrix [" << i << "] : " << endl
         << endl;
    for (int i = 0; i < el.len(); i++)
    {
        for (int j = 0; j < el.len(); j++)
            cout << scientific << el.getMatrixValue(i, j) << "  ";
        cout << endl;
    }
    cout << endl;
    cout << "Load Vector [" << i << "] : " << endl
         << endl;
    for (int i = 0; i < el.len(); i++)
        cout << scientific << el.getLoadVector(i) << "  ";
    cout << endl;
}

void assemble(BElement1D **el, matGlbM, double *GlbV, int n, int nElem)
{
    // assuming the Global matrix and global load vector are alloc'd
    // first we make sure they are set to 0
    for (int i = 0; i < n * nElem + 1; i++)
    {
        for (int j = 0; j < n * nElem + 1; j++)
            GlbM[i][j] = 0.0;
        GlbV[i] = 0.0;
    }
    // then we "assemble"
    for (int k = 0; k < nElem; k++)
    {
        for (int i = 0; i <= n; i++)
        {
            for (int j = 0; j <= n; j++)
                GlbM[n * k + i][n * k + j] += el[k]->getMatrixValue(i, j);

            GlbV[n * k + i] += el[k]->getLoadVector(i);
        }
    }
} // O(n^2 * nElem)

// difference by maximum norm
double vec_dif(double *v1, double *v2, int len)
{
    double max = -1.0e66;
    for (int i = 0; i < len; i++)
    {
        double dif = abs(v2[i] - v1[i]);
        if (dif > max)
            max = dif;
    }
    return max;
}
/*
// Jacobi iterative linear solver
double *linearSolve(matA, double *b, int len)
{
    double *res_k1 = new double[len];
    double *res_k2 = new double[len];
    double tolerance = 0.1e-5;
    int k;

    // sets res_k2 to b
    for (int i = 0; i < len; i++)
        res_k1[i] = b[i];

    // Gauss-Seidel iteration
    for (k = 0; k < MAX_ITERATIONS && vec_dif(res_k1, res_k2, len) > tolerance; k++)
    {
        // copy res_k to res_(k-1) (res_k2 to res_k1)
        for (int i = 0; i < len; i++)
            res_k1[i] = res_k2[i];

        // then compute res_k
        for (int i = 0; i < len; i++)
        {
            double s = 0.0;
            for (int j = 0; j < i; j++)
            {
                s += A[i][j] * res_k1[j];
            }
            for (int j = i+1; j < len; j++)
            {
                s += A[i][j] * res_k1[j];
            }
            res_k2[i] = (b[i] - s) / A[i][i];
        }
    }

    if (k == MAX_ITERATIONS)
        cerr << "Linear solve max iterations reached!" << endl
             << "Solution might show wrong results!" << endl;

    return res_k2;
}
*/

// Gaussian elimination linear solver
double *linearSolve(matA, double *b, int len)
{
    double *res = new double[len];
    
    // Makes the elimination
    for (int k = 0; k < len-1; k++)
    {
        for (int i = k+1; i < len; i++)
        {
            double m = A[i][k] / A[k][k];
            b[i] -= m * b[k];
            for (int j = k+1; j < len; j++)
                A[i][j] -= m * A[k][j];
        }
    }

    // Solves the upper triangular system resultant from the routine above
    for (int i = len-1; i >= 0; i--)
    {
        double s = 0.0;
        for (int j = i + 1; j < len; j++)
            s += res[j] * A[i][j];
        res[i] = (b[i] - s) / A[i][i];
    }

    return res;
}
