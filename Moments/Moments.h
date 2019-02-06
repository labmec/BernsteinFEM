#pragma once
#include <armadillo>
#include "Elements.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))

// Abstract class
template <typename _signature, Element_t EL>
class BMoment
{
  protected:
    int q;                // number of quadrature points in one dimension
    int n;                // Bernstein polynomial order
    int lenMoments;       // length of the Bmoment vector
    int lenCval;          // length of the Cval array
    arma::mat Bmoment;    // vector where the b-moments are stored
    arma::mat Cval;       // vector where the function values are stored
    bool fValSet = false; // is true if the function value is set
    bool fDefSet = false; // is true if the function definition is set
    bool functVal = true; // determines if will use function value or function definition (use function value by default)
    int nb_Array;         // dimension of function Image (function is scalar valued by default)
    Element<EL> element;  // element to be used (either 1D, quadrilater, triangle, etc.)
    arma::mat quadraWN;   // quadrature points and weights

    // function definition for the computation of the b-moments
    std::function<_signature> f;

    // protected virtual methods
    // assign the quadrature points and weights
    virtual void assignQuadra() = 0;

    // loads the function definition values at quadrature points into Cval
    virtual void loadFunctionDef() = 0;

  public:
    // constructors
    BMoment(int q, int n, const Element<EL> &element = Element<EL>(), int nb_Array = 1);

    // copy constructor
    BMoment(const BMoment<_signature, EL> &cp);

    // copy assignment operator
    BMoment &operator=(const BMoment &cp);

    // destructor
    ~BMoment();

    // getters
    // returns number of integration points
    int getNumIntegrationPoints();

    // returns polynomial order
    int getPOrder();

    // returns moments array length
    int getLenMoments();

    // returns number
    int getLenCval();

    // return the dimension of function Image (function is scalar valued by default)
    int getNbArray();

    // returns wether you will be using function values or function definition (true for function values)
    bool getFunctVal();

    // returns the element used for computing
    Element<EL> getElement();

    // returns the whole Bmoment matrix
    const arma::mat &getMoments();

    // setters
    // sets number of integration points
    void setNumIntegrationPoints(int q);

    // sets polynomial order
    void setPOder(int n);

    // sets the dimension of function Image (function is scalar valued by default)
    void setNbArray(int nb_Array);

    void setFunctionValues(const arma::vec &Cval);

    void setFunctionValues(const arma::mat &Cval);

    void setFunctionDefinition(std::function<_signature> f);

    void setElement(Element<EL> element);

    // inline methods

    void zero();

    void useFunctionDef();

    void useFunctionValue();

    void computeMoments(std::function<_signature> f);

    void computeMoments(const arma::mat &Cval);

    // returns the vector with the integration points over the object's element, following the moments organization
    virtual arma::mat getIntegrationPoints() = 0; // TODO: consider changing this to const arma::mat &

    // computes the moments and store it in the Bmoment array, use getBMoment to get it
    virtual arma::mat &computeMoments() = 0;
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
    BMoment1D(int q, int n, const Element<Element_t::LinearEl> &element = Element<Element_t::LinearEl>(), int nb_Array = 1);

    // copy constructor
    BMoment1D(const BMoment1D &cp);

    // copy assignment operator
    BMoment1D &operator=(const BMoment1D &cp);

    // destructor
    ~BMoment1D();

    // returns the value of the i-th indexed B-moment
    double getBMoment(int i) { return Bmoment(i, 0); }

    // returns the value of the i-th indexed B-moment at the specified dimension 'dim'
    double getBMoment(int i, int dim) { return Bmoment(i, dim - 1); }

    arma::mat getIntegrationPoints();

    // compute the B-moments using the values already assigned in the object
    arma::mat &computeMoments() final;
};

/*****************************************************************************
 * Bernstein moments for triangles/simpleces (2-dimensional)                 *
 *****************************************************************************/
class BMoment2DTri : public BMoment<double(double, double), Element_t::TriangularEl>
{
    arma::mat BMomentInter;

  protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    void loadFunctionDef() final;

  public:
    // default constructor
    BMoment2DTri(int q, int n, const Element<Element_t::TriangularEl> &element = Element<Element_t::TriangularEl>(), int nb_Array = 1);

    // copy constructor
    BMoment2DTri(const BMoment2DTri &cp);

    // copy assignment operator
    BMoment2DTri &operator=(const BMoment2DTri &cp);

    ~BMoment2DTri();

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1) on the specified dimension
    double getBMoment(uint a1, uint a2, int dim)
    {
        try
        {
            return Bmoment(element.position({a1, a2}, n), dim);
        }
        catch (std::logic_error &e)
        {
            std::cerr << "at function: " << __func__ << std::endl;
            std::cerr << "arguments values ot of range: " << std::endl;
            std::cerr << "a1 = " << a1 << "; a2 = " << a2 << "; dim = " << dim << std::endl;
            std::cerr << "a1 max: " << n << "a2 max:" << n - a1 << "dim max: " << nb_Array - 1 << std::endl;
            throw std::logic_error(e);
        }
    }

    // get the i-th bmoment in the array, only use if you really know what you're doing
    double getBMoment(int i) { return Bmoment(i, 0); }

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1)
    double getBMoment(int i, int dim) { return Bmoment(i, dim); }

    // returns the vectors with the integration points (x, y) over the object's element, following the moments organization
    // Assuming: points = getIntegrationPoints(); then
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    arma::mat getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    arma::mat &computeMoments() final;
};

/*****************************************************************************
 * Bernstein moments for quadrilaterals (2-dimensional)                      *
 *****************************************************************************/

// for later implementation:
// add the second order of polynomial degree to computation
class BMoment2DQuad : public BMoment<double(double, double), Element_t::QuadrilateralEl>
{
    int m;                  // second polynomial order
    arma::mat BMomentInter; // auxiliary matrix to compute moments

  protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    // computes the function by the definition and stores it in Cval
    void loadFunctionDef() final;

  public:
    // default constructor
    BMoment2DQuad(int q, int n, const Element<Element_t::QuadrilateralEl> &element = Element<Element_t::QuadrilateralEl>(), int nb_Array = 1);

    BMoment2DQuad(int q, int n, int m, const Element<Element_t::QuadrilateralEl> &element = Element<Element_t::QuadrilateralEl>(), int nb_Array = 1);

    // copy constructor
    BMoment2DQuad(const BMoment2DQuad &cp);

    // copy assignment operator
    BMoment2DQuad &operator=(const BMoment2DQuad &cp);

    // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral on the specified dimension
    double getBMoment(int i, int dim) { return Bmoment(i, dim); }

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 on the specified dimension
    double getBMoment(uint a1, uint a2, int dim) { return Bmoment(element.position({a1, a2}, n), dim); }

    // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral
    double getBMoment(int i) { return Bmoment(i, 0); }

    // returns the vector with the integration points (x, y) over the object's element, following the moments organization
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    arma::mat getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    arma::mat &computeMoments() final;
};

// Future class definition
class BMoment3DCube : public BMoment<double(double, double, double), Element_t::CubeEl>
{
    int m;                  // second polynomial order
    int p;                  // third polynomial order
    arma::mat BMomentInter; // auxiliary matrix to compute moments

  protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    // computes the function by the definition and stores it in Cval
    void loadFunctionDef() final;

  public:
    // default constructor
    BMoment3DCube(int q, int n, const Element<Element_t::CubeEl> &element = Element<Element_t::CubeEl>(), int nb_Array = 1);

    BMoment3DCube(int q, int n, int m, int p, const Element<Element_t::CubeEl> &element = Element<Element_t::CubeEl>(), int nb_Array = 1);

    // copy constructor
    BMoment3DCube(const BMoment3DCube &cp);

    // copy assignment operator
    BMoment3DCube &operator=(const BMoment3DCube &cp);

    // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral on the specified dimension
    double getBMoment(int i, int dim) { return Bmoment(i, dim); }

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 on the specified dimension
    double getBMoment(uint a1, uint a2, uint a3, int dim) { return Bmoment(element.position({a1, a2, a3}, n), dim); }

    // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral
    double getBMoment(int i) { return Bmoment(i, 0); }

    // returns the vector with the integration points (x, y) over the object's element, following the moments organization
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    arma::mat getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    arma::mat &computeMoments() final;
};

class BMoment3DTetra : public BMoment<double(double, double, double), Element_t::TetrahedronEl>
{
    arma::mat BMomentInter;

  protected:
    // map to obtain Gauss-Jacobi rule on unit interval
    void assignQuadra() final;

    void loadFunctionDef() final;

  public:
    // default constructor
    BMoment3DTetra(int q, int n, const Element<Element_t::TetrahedronEl> &element = Element<Element_t::TetrahedronEl>(), int nb_Array = 1);

    // copy constructor
    BMoment3DTetra(const BMoment3DTetra &cp);

    // copy assignment operator
    BMoment3DTetra &operator=(const BMoment3DTetra &cp);

    ~BMoment3DTetra();

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1) on the specified dimension
    double getBMoment(uint a1, uint a2, uint a3, int dim)
    {
        try
        {
            return Bmoment(element.position({a1, a2, a3}, n), dim);
        }
        catch (std::logic_error &e)
        {
            std::cerr << "at function: " << __func__ << std::endl;
            std::cerr << "arguments values ot of range: " << std::endl;
            std::cerr << "a1 = " << a1 << "; a2 = " << a2 << "; dim = " << dim << std::endl;
            std::cerr << "a1 max: " << n << "a2 max:" << n - a1 << "dim max: " << nb_Array - 1 << std::endl;
            throw std::logic_error(e);
        }
    }

    // get the i-th bmoment in the array, only use if you really know what you're doing
    double getBMoment(int i) { return Bmoment(i, 0); }

    // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1)
    double getBMoment(int i, int dim) { return Bmoment(i, dim); }

    // returns the vectors with the integration points (x, y) over the object's element, following the moments organization
    // Assuming: points = getIntegrationPoints(); then
    // points(i, 0) == x i-th coordinate
    // points(i, 1) == y i-th coordinate
    arma::mat getIntegrationPoints();

    // compute the b-moments using the values already assigned in the object
    arma::mat &computeMoments() final;
};