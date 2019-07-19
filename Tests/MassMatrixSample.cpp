/** MassMatrixSample
 * Author: Lucas Barretto Andrade
 *
 * This file contains a sample code showing how to compute
 * an element's Mass Matrix in Finite Element Analysis using
 * this library and the Bernstein polynomials
 */

#include <iostream>

// uncomment the next line to increase performance a little by disabling bounds check on array's element access
//#define ARMA_NO_DEBUG
#include "BernsteinFEM.h"

using namespace BernsteinFEM;

double function(double);

int main()
{
	int polynomialOrder = 2;	// basis polynomial order
	int numIntegrationPoints = 2 * polynomialOrder;	// number of integration points

	TPZFMatrix<REAL> vertices({ 0.0, 1.0 }); // vertices of the element (default: {0.0, 1.0})

	TPZVec<uint64_t> indexVector({ 0, 1 }); // vector with the indices of the element's vertices (default {0, 1})

	Element<Element_t::LinearEl> element(vertices, indexVector); // element in question

	BMass1D MassMatrixComputer(numIntegrationPoints, polynomialOrder, element); // object that manages Mass Matrix computation

	auto integrationPoints = MassMatrixComputer.getIntegrationPoints(); // gets the vector with integration points

	TPZVec<REAL> functionValues(integrationPoints.Rows());
	
	// evaluate the function in integration points and store it in 'functionValues'
	for (unsigned i = 0; i < functionValues.size(); i++)
	{
		functionValues[i] = function(integrationPoints(i));
	}

	// pass the values to the MassMatrixComputer
	MassMatrixComputer.setFunctionValues(functionValues);
	// additionaly you can just use the functor, like this:
	// MassMatrixCompute.setFunctionDefinition(function)

	// compute the Mass Matrix
	MassMatrixComputer.computeMatrix();

	// gets the mass matrix object in a actual matrix (arma::mat)
	auto massMatrix = MassMatrixComputer.getMatrix();

	// from this point on you can delete the BMass1D object if managing memory dinamically

	massMatrix.Print("Mass Matrix:", std::cout, EFormatted);
}

// sample function
double function(double x)
{
	return 1.0;
}