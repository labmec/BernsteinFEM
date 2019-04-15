/** mom-tri.cpp
 * 
 * This file makes tests for the computation
 * of Moments of an 2D simplex element using
 * Bernstein polynomials.
 * 
 * A Moment is the integral of a basis function
 * (Bernstein polynomials) multiplied by a weight function.
 */

#include <iostream>

#define _USE_MATH_DEFINES
#include <cmath>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include <armadillo>
#include "Moments.h"
#include "JacobiGaussNodes.h"

using namespace std;

// testing if BMoments is functioning correctly
// you can change the function to use in the computations by
// modifying the function below, and uncommenting the code under
// the "sets function values" section

double f(double x)
{
    return sin(M_PI * x);
}

int main()
{
    int n = 4;
    int q = 2 * n;

    BMoment2DTri mom_tri(q, n);

    // sets triangle to compute
    Element<Element_t::TriangularEl> tri_el({{0,0},{0,1},{1,0}});
    mom_tri.setElement(tri_el);

    // this last part is only necessary when computing from function definition

    // sets function values
    {
        arma::vec Fval((q + 1) * (q + 1), arma::fill::ones); // it has ((q + 2) * (q + 1)) / 2 nodes, but is indexed using (q + 1) * (q + 1) positions

        // for (int i = 0; i < q * q; i++)
        //     Fval(i) = f( (legendre_xi(q, i) + 1.0) * 0.5 );

        mom_tri.setFunctionValues(Fval);
    }

    // computes moments
    mom_tri.computeMoments();

    // prints moments
    cout << "Bernstein moments computed with the function sin(pi * x)" << endl
         << endl;
    try 
    {
        for (int a1 = 0; a1 < (n + 1); a1++)
        for (int a2 = 0; a2 <= n - a1; a2++)
            cout << "Bmoment[" << a1 << ", " << a2 << ", " << n - a1 - a2 << "] : " 
            << mom_tri.getBMoment(a1, a2, 0) /* this takes indexing into account */
            << endl;
    }
    catch(std::logic_error &e)
    {
        cerr << "at file: "<< __FILE__ << ", function: " << __func__ << endl;
    }

    return 0;
}