/** LoadVecSample.cpp
 * Author: Lucas Barretto Andrade
 * 
 * This file contains a sample code showing how to compute
 * an element's Load Vector in Finite Element Analysis using
 * this library and the Bernstein polynomials
 */

#include <iostream>

// uncomment the next line to increase performance a little by disabling bounds check on array's element access
//#define ARMA_NO_DEBUG
#include "BernsteinFEM.h"

using namespace BernsteinFEM;

// load function
double function(double x);

int main ()
{
    int polynomialOrder = 2;        // basis polynomial order
    int nIntegrationPoints = 2 * polynomialOrder; // number of integration points

    arma::vec vertices({0.0, 1.0}); // vertices of the element in question (default = {0.0, 1.0})

    arma::ivec indexVector({0, 1}); // vector with the indices of the element's vertices (default = {0, 1})

    Element<Element_t::LinearEl> element(vertices, indexVector); // element in question, built with vertices and its indices
    // omitting the arguments in the constructor above would give it the default values, shown above

    BLoad1D LoadVecComputer(nIntegrationPoints, polynomialOrder, element); // object that manages a load vector computation

    auto integrationPoints = LoadVecComputer.getIntegrationPoints(); // gets the vector with integration points for you to evaluate your load function into

    arma::vec functionEvaluated(integrationPoints.n_rows); // creates the vector which will hold the values for the evaluation of the function in the integration points

    // evaluates the function in the integration points and store it in 'functionEvaluated'
    for(unsigned i = 0; i < integrationPoints.n_rows; i++)
    {
        functionEvaluated(i) = function(integrationPoints(i));
    }

    // pass the values to the LoadVecComputer object
    LoadVecComputer.setFunctionValues(functionEvaluated);

    // compute the Load Vector
    auto loadVector = LoadVecComputer.computeMoments();

	// from this point on you can delete the BLoad1D object if managing memory dinamically

    // prints the vector using the built-in armadillo output
    loadVector.print("Load Vector");
}

// sample load function
double function(double x)
{
    return x;
}