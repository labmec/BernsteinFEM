#pragma once
#include <armadillo>
#include "Elements.h"

#define MAX(a, b) ((a) < (b) ? (b) : (a))

// Abstract class
template <typename _signature, Element_t EL>
class BMoment
{
  int q; // number of quadrature points in one dimension
  int n; // Bernstein polynomial order

protected:
  int lenMoments;       // length of the Bmoment vector
  int lenCval;          // length of the Cval array
  arma::mat Bmoment;    // vector where the b-moments are stored
  arma::mat Cval;       // vector where the function values are stored
  bool fValSet = false; // is true if the function value is set
  bool fDefSet = false; // is true if the function definition is set
  int functVal;         // determines if will use function value or function definition (use function dvalue by default)
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
  BMoment(int q, int n, int nb_Array = 1, Element<EL> element = Element<EL>())
    : element(element), Bmoment(), Cval(), quadraWN()
  {
    this->q = q;
    this->n = n;
    this->nb_Array = nb_Array;
  }

  // copy constructor
  BMoment(const BMoment<double(double), Element_t::LinearEl> &cp)
  {
    this->q = cp.q;
    this->n = cp.n;
    this->nb_Array = cp.nb_Array;
  }

// destructor
  ~BMoment();

// getters
  // returns number of integration points
  int getNumIntegrationPoints() { return q; }

  // returns polynomial order
  int getPOrder() { return n; }

  // returns moments array length
  int getLenMoments() { return lenMoments; }

  // returns number
  int getLenCval() { return lenCval; }

  // return the dimension of function Image (function is scalar valued by default)
  int getNbArray() { return nb_Array; }

  // returns wether you will be using function values or function definition (0 or 1, respectively)
  int getFunctVal() { return functVal; }

  // returns the element used for computing
  Element<EL> getElement() { return element; }

  // returns the whole Bmoment matrix
  const arma::mat &getBMoment() { return Bmoment; }

// setters
  // sets number of integration points
  void setNumIntegrationPoints(int q) { this->q = q; }

  // sets polynomial order
  void setPOder(int n) { this->n = n; }

  // sets the dimension of function Image (function is scalar valued by default)
  void setNbArray(int nb_Array)
  {
    this->nb_Array = nb_Array;
    Bmoment.resize(lenMoments, nb_Array);
    Cval.resize(lenCval, nb_Array);
  }

  void setFunctionValues(const arma::vec &Cval) { this->Cval.swap(Cval); }

  void setFunctionValues(const arma::mat &Cval) { this->Cval.swap(Cval); }

  void setFunctionDefinition(std::function<_signature> f) { this->f = f; }

  void setElement(Element<EL> element) { this->element = element; }

// inline methods

  void zero() { Bmoment.zeros(); }

  void useFunctionDef() { functVal = 0; }

  void useFunctionValue() { functVal = 1; }

  void compute_moments(std::function<_signature> f)
  {
    setFunctionDefinition(f);
    compute_moments();
  }

  void compute_moments(const arma::mat &Cval)
  {
    setFunctionValues(Cval);
    compute_moments();
  }

// virtual methods
  // returns the position of a specified point in the matrix (the point is defined by the element)
  virtual int position() = 0;

  // returns the vector with the integration points over the object's element, following the moments organization
  virtual arma::vec getIntegrationPoints() = 0;

  // computes the moments and store it in the Bmoment array, use getBMoment to get it
  virtual void computeMoments() = 0;
};

/*****************************************************************************
 * 1-dimensional Bernstein moments                                           *
 *****************************************************************************/
class BMoment1D : BMoment<double(double), Element_t::LinearEl>
{
public:
  // constructors
  BMoment1D(int q, int n, int nb_Array = 1, Element<Element_t::LinearEl> element = Element<Element_t::LinearEl>())
    : BMoment(q, n, nb_Array, element) { }

  // destructor
  ~BMoment1D();

  // returns the index of the i-th index on the interval (unnecessary in this case, just made to be concise with the other versions)
  int position(int i, int n) { return i; }

  // zeroes the moments vector
  void zero() { Bmoment.zeros(); }

  // returns the whole Bmoment matrix
  const arma::mat &get_bmoment() { return Bmoment; }

  // returns the value of the i-th indexed B-moment
  double get_bmoment(int i) { return Bmoment(i, 0); }

  // returns the value of the i-th indexed B-moment at the specified dimension 'dim'
  double get_bmoment(int i, int dim) { return Bmoment(i, dim - 1); }

  // returns the vector with the integration points over the object's element, following the moments organization
  arma::vec getIntegrationPoints();

  // call if you're going to use the function definition as parameters instead of the function value (as in default)
  void useFunctionDef() { functVal = 0; }

  // call if you're going back to using the function values
  void useFunctionValue() { functVal = 1; }

  // set the function values for computation
  void setFunction(const arma::vec &Fval);

  // set the function values for computation when you have more than one dimension
  void setFunction(const arma::mat &Fval);

  // set the function definition for computation
  void setFunction(std::function<double(double)> f);

  void setInterval(double a, double b);

  // compute the B-moments using the values already assigned in the object
  void compute_moments();

  // compute the b-moments for the specified f function
  void compute_moments(std::function<double(double)> f);

  // compute the b-moments for the Fval function values
  void compute_moments(const arma::vec &Fval);
};

/*****************************************************************************
 * Bernstein moments for triangles/simpleces (2-dimensional)                 *
 *****************************************************************************/
class BMoment2DTri
{
  int q;                // number of quadrature points in one dimension
  int n;                // Bernstein polynomial order
  int lenMoments;       //length of the Bmoment vectors
  arma::mat CVal;       // stores function values at quadrature points
  arma::mat Bmoment;    // Vectors to store the Bernstein Moments
  bool fValSet = false; // is true if the function value is set
  bool fDefSet = false; // is true if the function definition is set

  // function definition for the computation of the b-moments
  std::function<double(double, double)> f;

  // map to obtain Gauss-Jacobi rule on unit interval
  void assignQuadra();

  //convert barycentric coordinates (b1,b2,b3) of a point w.r.t. vertices v1,v2,v3 into Cartesian coordinates v
  void bary2cart2d(double b1, double b2, double b3, double v1[2], double v2[2], double v3[2], double v[2]);

  //initialize Bmoment by the values of the function f at the quadrature points of order q
  void init_BmomentC_Bmom2d();

  //initialize Bmoment with the values of the function f, stored at CVal
  void init_Bmoment2d_Cval();

protected:
  int functVal = 1;   // determines if will use function value or function definition (use function value by default)
  int nb_Array = 1;   // dimension of function Image (function is scalar valued by default)
  arma::mat vertices; // triangle vertices coordinates
  arma::mat quadraWN; // quadrature points and weights

  // helps indexing quadrature points vectors
  int position_q(int i, int j, int q) { return i * q + j; }

  // routine used to pre-multiply normals with the Bmoments
  void transform_BmomentC_Stiff2d(BMoment2DTri *Bmomentab, const arma::mat &normalMat);

public:
  // default constructor
  BMoment2DTri();

  // quadrature and polynomial order constructor;
  BMoment2DTri(int q, int n);

  BMoment2DTri(int q, int n, int nb_Array);

  // constructor setting q, n and the triangle vertices
  BMoment2DTri(int q, int n, double T[][2]);

  BMoment2DTri(int q, int n, int nb_Array, double T[][2]);

  ~BMoment2DTri();

  // computes area of triangle < v1,v2,v3 >
  static double Area2d(double v1[2], double v2[2], double v3[2])
  {
    double x1 = v1[0];
    double y1 = v1[1];
    double x2 = v2[0];
    double y2 = v2[1];
    double x3 = v3[0];
    double y3 = v3[1];

    return abs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2;
  }

  // computes area of triangle defined by vertices
  static double Area2d(const arma::mat &vertices)
  {
    double x1 = vertices(0, 0);
    double y1 = vertices(0, 1);
    double x2 = vertices(1, 0);
    double y2 = vertices(1, 1);
    double x3 = vertices(2, 0);
    double y3 = vertices(2, 1);

    return fabs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2.;
  }

  // sets the dimension number of the function multiplying the Bernstein basis polynomial
  void setNbArray(int nb_Array)
  {
    this->nb_Array = nb_Array;
    Bmoment.resize(lenMoments, nb_Array);
    CVal.resize(MAX(n + 1, q) * MAX(n + 1, q), nb_Array);
  }

  int getNbArray() { return nb_Array; }

  // return the index for the (i, j, n-i-j) triangle coordinate
  static int position(int i, int j, int n) { return i * (n + 1) + j; }

  // zeroes the moments vector
  void zero() { Bmoment.zeros(); }

  // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1) on the specified dimension
  double get_bmoment(int a1, int a2, int dim) { return Bmoment(position(a1, a2, n), dim); }

  // get the i-th bmoment in the array, only use if you really know what you're doing
  double get_bmoment(int i) { return Bmoment(i, 0); }

  // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1)
  double get_bmoment(int i, int dim) { return Bmoment(i, dim); }

  // returns the vectors with the integration points (x, y) over the object's element, following the moments organization
  // Assuming: points = getIntegrationPoints(); then
  // points(i, 0) == x i-th coordinate
  // points(i, 1) == y i-th coordinate
  arma::mat getIntegrationPoints();

  // call if you're going to use the function definition as parameters instead of the function value (as in default)
  void useFunctionDef() { functVal = 0; }

  // call if you're going back to using the function values
  void useFunctionValue() { functVal = 1; }

  // set the function value at quadrature points, as in Fval, the Fval vector must use the order given by the position() function
  void setFunction(const arma::vec &Fval);
  void setFunction(const arma::mat &Fval);

  // set the function that multiplies the B-polynomial by definition
  void setFunction(std::function<double(double, double)> f);

  // set the element triangle vertices
  void setTriangle(double v1[2], double v2[2], double v3[2]);
  void setTriangle(const arma::mat &vertices);

  // computes the function definition into the function values vector
  void computeFunctionDef();

  // compute the b-moments using the values already assigned in the object
  void compute_moments();

  // compute the b-moments for the specified f function
  void compute_moments(std::function<double(double, double)> f);

  // compute the b-moments for the Fval function values
  void compute_moments(const arma::vec &Fval);
};

/*****************************************************************************
 * Bernstein moments for quadrilaterals (2-dimensional)                      *
 *****************************************************************************/

// for later implementation:
// add the second order of polynomial degree to computation
class BMoment2DQuad
{
  int q;                // number of quadrature points in one dimension
  int n;                // Bernstein polynomial order
  int m;                // Bernstein polynomial order (second independent variable)
  int lenMoments;       // length of the Bmoment vector
  arma::mat Bmoment;    // vector where the b-moments are stored
  arma::mat Cval;       // vector where the function values are stored
  bool fValSet = false; // is true if the function value is set
  bool fDefSet = false; // is true if the function definition is set

  // function definition for the computation of the b-moments
  std::function<double(double, double)> f;

  // methods

  // map to obtain Gauss-Jacobi rule on unit interval
  void assignQuadra();

  // computes the function by the definition and stores it in Cval
  void computeFunctionDef();

protected:
  int functVal = 1;   // determines if will use function value or function definition (use function value by default)
  int nb_Array = 1;   // dimension of function Image (function is scalar valued by default)
  arma::mat quadraWN; // quadrature points and weights
  arma::mat vertices; // vertices of the quadrilateral element

  // helps indexing quadrature points vectors
  int position_q(int i, int j, int q) { return i * q + j; }

public:
  // constructors
  BMoment2DQuad();
  BMoment2DQuad(int q, int n);
  BMoment2DQuad(int q, int n, int m, int nb_Array = 1);

  // destructor
  ~BMoment2DQuad();

  // sets the dimension number of the function multiplying the Bernstein basis polynomial
  void setNbArray(int nb_Array)
  {
    this->nb_Array = nb_Array;
    Bmoment.resize(lenMoments, nb_Array);
    Cval.resize(q * q, nb_Array);
  }

  int getNbArray() { return nb_Array; }

  // return the index for the (i, j) quadrilateral node
  static int position(int i, int j, int n) { return i * (n + 1) + j; }

  // zeroes the moments vector
  void zero() { Bmoment.zeros(); }

  // return the Bmoment matrix
  const arma::mat &get_bmoment() { return Bmoment; }

  // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral on the specified dimension
  double get_bmoment(int i, int dim) { return Bmoment(i, dim); }

  // get the bmoment value of the Bernstein polynomial with indexes a1 and a2 on the specified dimension
  double get_bmoment(int a1, int a2, int dim) { return Bmoment(position(a1, a2, n), dim); }

  // get the i-th Bmoment in the array, associated with the i-th node of the quadrilateral
  double get_bmoment(int i) { return Bmoment(i, 0); }

  // returns the vector with the integration points (x, y) over the object's element, following the moments organization
  // points(i, 0) == x i-th coordinate
  // points(i, 1) == y i-th coordinate
  arma::mat getIntegrationPoints();

  // call if you're going to use the function definition as parameters instead of the function value (as in default)
  void useFunctionDef() { functVal = 0; }

  // call if you're going back to using the function values
  void useFunctionValue() { functVal = 1; }

  // set the function value at quadrature points, as in Fval, the Fval vector must use the order given by the position() function
  void setFunction(const arma::vec &Fval);
  void setFunction(const arma::mat &Fval);

  // set the function that multiplies the B-polynomial
  void setFunction(std::function<double(double, double)> f);

  // set the element quadrilateral vertices
  void setQuadrilateral(double v1[2], double v2[2], double v3[2], double v4[2]);
  void setQuadrilateral(const arma::vec &v1, const arma::vec &v2, const arma::vec &v3, const arma::vec &v4);
  void setQuadrilateral(const arma::mat &vertices);

  // Returns the nodal shape function of the elements quadrilateral, and the jacobian determinant
  // the points are stored in the parameter X and the jacobian determinant in dX
  void nodalShape(double X[], double &dX, double xi, double eta);
  void nodalShape(arma::vec &X, double &dX, double xi, double eta);
  void nodalShape(arma::vec &X, arma::mat &jac, double &dX, double xi, double eta);

  // compute the b-moments using the values already assigned in the object
  void compute_moments();

  // compute the b-moments for the specified f function
  void compute_moments(std::function<double(double, double)> f);

  // compute the b-moments for the Fval function values
  void compute_moments(const arma::vec &Fval);
};

// Future class definition
class BMoment3D
{
};