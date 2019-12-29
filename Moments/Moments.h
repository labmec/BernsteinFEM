#pragma once

#include "pzreal.h"	  // REAL, STATE
#include "pzvec.h"	  // Vec
#include "pzmatrix.h" // TPZFMatrix
#include "Elements.h"
#include <vector>
#include <functional>

#define MAX(a, b) ((a) < (b) ? (b) : (a))

// auxiliary functions
// computes area of triangle < v1,v2,v3 >
inline
double triangle_area(double v1[2], double v2[2], double v3[2])
{
	double x1 = v1[0];
	double y1 = v1[1];
	double x2 = v2[0];
	double y2 = v2[1];
	double x3 = v3[0];
	double y3 = v3[1];

	return fabs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2.;
}

// computes area of triangle defined by vertices
inline
double triangle_area(TPZFMatrix<REAL>& vertices)
{
	double x1 = vertices(0, 0);
	double y1 = vertices(0, 1);
	double x2 = vertices(1, 0);
	double y2 = vertices(1, 1);
	double x3 = vertices(2, 0);
	double y3 = vertices(2, 1);

	return fabs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2.;
}

// BMoment Interface
class BMomentT
{
protected:
	// protected virtual methods
	// assign the quadrature points and weights
	virtual void assignQuadra() = 0;

	// loads the function definition values at quadrature points into Cval
	virtual void loadFunctionDef() = 0;

public:
	// getters
	// returns number of integration points
	virtual uint getNumIntegrationPoints() = 0;

	// returns polynomial order
	virtual uint getPOrder() = 0;

	// returns moments array length
	virtual uint getLenMoments() = 0;

	// returns number
	virtual uint getLenCval() = 0;

	// returns wether you will be using function values or function definition (true for function values)
	virtual bool getFunctVal() = 0;

	// returns the whole Bmoment matrix
	virtual const TPZVec<REAL>& getMoments() = 0;

	// setters
	// sets number of integration points
	virtual void setNumIntegrationPoints(uint q) = 0;

	// sets polynomial order
	virtual void setPOder(uint n) = 0;

    // copies values of the load function evaluated at the integration points into this object
	virtual void setFunctionValues(const TPZVec<REAL>& Cval) = 0;

    // moves vector pointer instead of copying the values
    virtual void setFunctionValues(TPZVec<REAL> &&Cval) = 0;

	// inline methods

	virtual void zero() = 0;

	virtual void useFunctionDef() = 0;

	virtual void useFunctionValue() = 0;

	virtual void computeMoments(const TPZVec<REAL>& Cval) = 0;

	// returns the vector with the integration points over the object's element, following the moments organization
	virtual TPZFMatrix<REAL> getIntegrationPoints() = 0; // TODO: consider changing this to const TPZVec<REAL> &

	// computes the moments and store it in the Bmoment array, use getBMoment to get it
	virtual TPZVec<REAL>& computeMoments() = 0;
};

// BMoment Abstract templated class
template <typename _signature, Element_t EL>
class BMoment : public BMomentT
{
protected:
    uint q;                  // number of quadrature points in one dimension
    uint n;                  // Bernstein polynomial order
    uint lenMoments;         // length of the Bmoment vector
    uint lenCval;            // length of the Cval array
    TPZVec<REAL> Bmoment;	 // vector where the b-moments are stored
    TPZVec<REAL> Cval;       // vector where the function values are stored
    bool fValSet = false;    // is true if the function value is set
    bool fDefSet = false;    // is true if the function definition is set
    bool functVal = true;    // determines if will use function value or function definition (use function value by default)
    Element<EL> element;     // element to be used (either 1D, quadrilater, triangle, etc.)
    TPZVec<REAL> intPoints;  // quadrature points
    TPZVec<REAL> intWeights; // quadrature weights

    // function definition for the computation of the b-moments
    std::function<_signature> f;

    // protected virtual methods
    // assign the quadrature points and weights
    virtual void assignQuadra() = 0;

    // loads the function definition values at quadrature points into Cval
    virtual void loadFunctionDef() = 0;

public:
    // constructors
    BMoment(uint q, uint n, const Element<EL> &element = Element<EL>());

    // copy constructor
    BMoment(const BMoment<_signature, EL> &cp);

    // copy assignment operator
    BMoment &operator=(const BMoment &cp);

    // destructor
    virtual ~BMoment() = default;

    // getters
    // returns number of integration points
    uint getNumIntegrationPoints();

    // returns polynomial order
    uint getPOrder();

    // returns moments array length
    uint getLenMoments();

    // returns number
    uint getLenCval();

    // returns wether you will be using function values or function definition (true for function values)
    bool getFunctVal();

    // returns the element used for computing
    Element<EL> &getElement();

    // returns the whole Bmoment matrix
    const TPZVec<REAL> &getMoments();

    // setters
    // sets number of integration points
    void setNumIntegrationPoints(uint q);

    // sets polynomial order
    void setPOder(uint n);

    void setFunctionValues(const TPZVec<REAL> &Cval);

    void setFunctionValues(TPZVec<REAL> &&Cval);

    void setFunctionDefinition(std::function<_signature> f);

    void setElement(Element<EL> element);

    // inline methods

    void zero();

    void useFunctionDef();

    void useFunctionValue();

    void computeMoments(std::function<_signature> f);

    void computeMoments(const TPZVec<REAL> &Cval);

    // returns the vector with the integration points over the object's element, following the moments organization
    virtual TPZFMatrix<REAL> getIntegrationPoints() = 0; // TODO: consider changing this to const TPZVec<REAL> &

    // computes the moments and store it in the Bmoment array, use getBMoment to get it
    virtual TPZVec<REAL> &computeMoments() = 0;
};

/*****************************************************************************
 * 1-dimensional Bernstein moments                                           *
 *****************************************************************************/
class BMoment1D : public BMoment<double(double), Element_t::LinearEl>
{
protected:
    void assignQuadra() final;

    void loadFunctionDef() final;

public:
    // constructors
    BMoment1D(uint q, uint n, const Element<Element_t::LinearEl> &element = Element<Element_t::LinearEl>());

    // copy constructor
    BMoment1D(const BMoment1D &cp);

    // copy assignment operator
    BMoment1D &operator=(const BMoment1D &cp);

    // destructor
    virtual ~BMoment1D() = default;

    // returns the value of the i-th indexed B-moment
    double getBMoment(uint i) { return Bmoment[i]; }

    TPZFMatrix<REAL> getIntegrationPoints();

    // compute the B-moments using the values already assigned in the object
    TPZVec<REAL> &computeMoments() final;
};

/*****************************************************************************
 * Bernstein moments for triangles/simpleces (2-dimensional)                 *
 *****************************************************************************/
class BMoment2DTri : public BMoment<double(double, double), Element_t::TriangularEl>
{
    TPZVec<REAL> BMomentInter;

protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    void loadFunctionDef() final;

public:
    // default constructor
    BMoment2DTri(uint q, uint n, const Element<Element_t::TriangularEl> &element = Element<Element_t::TriangularEl>());

    // copy constructor
    BMoment2DTri(const BMoment2DTri &cp);

    // copy assignment operator
    BMoment2DTri &operator=(const BMoment2DTri &cp);

    virtual ~BMoment2DTri() = default;

    // get the i-th bmoment in the array, only use if you really know what you're doing
    double getBMoment(uint i) { return Bmoment[i]; }

    // returns the vectors with the integration points (x, y) over the object's element, following the moments organization
    // Assuming: points = getIntegrationPoints(); then
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    TPZFMatrix<REAL> getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    TPZVec<REAL> &computeMoments() final;
};

/*****************************************************************************
 * Bernstein moments for quadrilaterals (2-dimensional)                      *
 *****************************************************************************/

// for later implementation:
// add the second order of polynomial degree to computation
class BMoment2DQuad : public BMoment<double(double, double), Element_t::QuadrilateralEl>
{
    uint m;                        // second polynomial order
    TPZVec<REAL> BMomentInter; // auxiliary matrix to compute moments

protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    // computes the function by the definition and stores it in Cval
    void loadFunctionDef() final;

public:
    // default constructor
    BMoment2DQuad(uint q, uint n, const Element<Element_t::QuadrilateralEl> &element = Element<Element_t::QuadrilateralEl>());

    BMoment2DQuad(uint q, uint n, uint m, const Element<Element_t::QuadrilateralEl> &element = Element<Element_t::QuadrilateralEl>());

    // copy constructor
    BMoment2DQuad(const BMoment2DQuad &cp);

	virtual ~BMoment2DQuad() = default;

    // copy assignment operator
    BMoment2DQuad &operator=(const BMoment2DQuad &cp);

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 on the specified dimension
    double getBMoment(uint a1, uint a2, int dim) { return Bmoment[element.position({a1, a2})]; }

    // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral
    double getBMoment(uint i) { return Bmoment[i]; }

    // returns the vector with the integration points (x, y) over the object's element, following the moments organization
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    TPZFMatrix<REAL> getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    TPZVec<REAL> &computeMoments() final;
};

// Future class definition
class BMoment3DCube : public BMoment<double(double, double, double), Element_t::CubeEl>
{
    uint m;                        // second polynomial order
    uint p;                        // third polynomial order
    TPZVec<REAL> BMomentInter; // auxiliary matrix to compute moments

protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    // computes the function by the definition and stores it in Cval
    void loadFunctionDef() final;

public:
    // default constructor
    BMoment3DCube(uint q, uint n, const Element<Element_t::CubeEl> &element = Element<Element_t::CubeEl>());

    BMoment3DCube(uint q, uint n, uint m, uint p, const Element<Element_t::CubeEl> &element = Element<Element_t::CubeEl>());

    // copy constructor
    BMoment3DCube(const BMoment3DCube &cp);

    // copy assignment operator
    BMoment3DCube &operator=(const BMoment3DCube &cp);

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 on the specified dimension
    double getBMoment(uint a1, uint a2, uint a3, int dim) { return Bmoment[element.position({a1, a2, a3})]; }

    // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral
    double getBMoment(uint i) { return Bmoment[i]; }

    // returns the vector with the integration points (x, y) over the object's element, following the moments organization
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    TPZFMatrix<REAL> getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    TPZVec<REAL> &computeMoments() final;
};

class BMoment3DTetra : public BMoment<double(double, double, double), Element_t::TetrahedronEl>
{
    TPZVec<REAL> BMomentInter;

protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    void loadFunctionDef() final;

public:
    // default constructor
    BMoment3DTetra(uint q, uint n, const Element<Element_t::TetrahedronEl> &element = Element<Element_t::TetrahedronEl>());

    // copy constructor
    BMoment3DTetra(const BMoment3DTetra &cp);

    // copy assignment operator
    BMoment3DTetra &operator=(const BMoment3DTetra &cp);

    ~BMoment3DTetra();

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1) on the specified dimension
    double getBMoment(uint a1, uint a2, uint a3, int dim)
    {
            return Bmoment[element.position({a1, a2, a3})];
    }

    // get the i-th bmoment in the array, only use if you really know what you're doing
    double getBMoment(uint i) { return Bmoment[i]; }

    // returns the vectors with the integration points (x, y) over the object's element, following the moments organization
    // Assuming: points = getIntegrationPoints(); then
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    TPZFMatrix<REAL> getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    TPZVec<REAL> &computeMoments() final;
};