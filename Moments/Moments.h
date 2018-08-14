#include <armadillo>

#ifndef MOMENTS_H
#define MOMENTS_H

#define MAX(a, b) ((a) < (b) ? (b) : (a))

/*****************************************************************************
 * 1-dimensional Bernstein moments                                           *
 *****************************************************************************/
class BMoment1D
{
  int q;                // number of quadrature points in one dimension
  int n;                // Bernstein polynomial order
  int lenMoments;       // length of the Bmoment vector
  arma::mat Bmoment;    // vector where the b-moments are stored
  arma::mat Cval;       // vector where the function values are stored
  bool fValSet = false; // is true if the function value is set
  bool fDefSet = false; // is true if the function definition is set

  // function definition for the computation of the b-moments
  double (*f)(double) = 0x0;

  // assign the quadrature points and weights
  void assignQuadra();

  // loads the function definition values at quadrature points into Cval
  void loadFunctionDef();

protected:
  int functVal = 1;   // determines if will use function value or function definition (use function dvalue by default)
  int nb_Array = 1;   // dimension of function Image (function is scalar valued by default)
  double a, b;        // interval [a, b] variables
  arma::mat quadraWN; // quadrature points and weights

public:
  // constructors
  // default constructor
  BMoment1D();
  // quadrature and polynomial order constructor
  BMoment1D(int q, int n);
  BMoment1D(int q, int n, int nb_Array);

  // destructor
  ~BMoment1D();

  // sets the dimension number of the function multiplying the Bernstein basis polynomial (default value: 1)
  void setNbArray(int nb_Array)
  {
    this->nb_Array = nb_Array;
    Bmoment.resize(lenMoments, nb_Array);
    Cval.resize(q * q, nb_Array);
  }

  int getNbArray() { return nb_Array; }

  // returns the index of the i-th index on the interval (unnecessary in this case, just made to be concise with the other versions)
  int position(int i, int n) { return i; }

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
  void setFunction(arma::vec Fval);
  
  // set the function values for computation when you have more than one dimension
  void setFunction(arma::mat Fval);

  // set the function definition for computation
  void setFunction(double (*f)(double));

  void setInterval(double a, double b);

  // compute the B-moments using the values already assigned in the object
  void compute_moments();

  // compute the b-moments for the specified f function
  void compute_moments(double (*f)(double));

  // compute the b-moments for the Fval function values
  void compute_moments(arma::vec Fval);
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
  double (*f)(double, double) = 0x0;

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

  // routine used to pre-multiply normals with the Bmoments
  void transform_BmomentC_Stiff2d(BMoment2DTri *Bmomentab, arma::mat normalMat);

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
  static double Area2d(arma::mat vertices)
  {
    double x1 = vertices(0, 0);
    double y1 = vertices(0, 1);
    double x2 = vertices(1, 0);
    double y2 = vertices(1, 1);
    double x3 = vertices(2, 0);
    double y3 = vertices(2, 1);

    return abs(x2 * y3 - x1 * y3 - x3 * y2 + x1 * y2 + x3 * y1 - x2 * y1) / 2;
  }

  // sets the dimension number of the function multiplying the Bernstein basis polynomial
  void setNbArray(int nb_Array)
  {
    this->nb_Array = nb_Array;
    Bmoment.resize(lenMoments, nb_Array);
    CVal.resize(q, nb_Array);
  }

  int getNbArray() { return nb_Array; }

  // return the index for the (i, j, n-i-j) triangle coordinate
  static int position(int i, int j, int n) { return i * (n + 1) + j; }

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
  void setFunction(arma::vec Fval);
  void setFunction(arma::mat Fval);

  // set the function that multiplies the B-polynomial by definition
  void setFunction(double (*function)(double, double));

  // set the element triangle vertices
  void setTriangle(double v1[2], double v2[2], double v3[2]);
  void setTriangle(arma::mat);

  // computes the function definition into the function values vector
  void computeFunctionDef();

  // compute the b-moments using the values already assigned in the object
  void compute_moments();

  // compute the b-moments for the specified f function
  void compute_moments(double (*f)(double, double));

  // compute the b-moments for the Fval function values
  void compute_moments(arma::vec Fval);
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
  int lenMoments;       // length of the Bmoment vector
  arma::mat Bmoment;    // vector where the b-moments are stored
  arma::mat Cval;       // vector where the function values are stored
  bool fValSet = false; // is true if the function value is set
  bool fDefSet = false; // is true if the function definition is set

  // function definition for the computation of the b-moments
  double (*f)(double, double) = 0x0;

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

public:
  // constructors
  BMoment2DQuad();
  BMoment2DQuad(int q, int n);
  BMoment2DQuad(int q, int n, int nb_Array);

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

  // return the Bmoment matrix
  const arma::mat& get_bmoment() { return Bmoment; }
  
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
  void setFunction(arma::vec Fval);
  void setFunction(arma::mat Fval);

  // set the function that multiplies the B-polynomial
  void setFunction(double (*function)(double, double));

  // set the element quadrilateral vertices
  void setQuadrilateral(double v1[2], double v2[2], double v3[2], double v4[2]);
  void setQuadrilateral(arma::vec v1, arma::vec v2, arma::vec v3, arma::vec v4);
  void setQuadrilateral(arma::mat vertices);

  // Returns the nodal shape function of the elements quadrilateral, and the jacobian determinant
  // the points are stored in the parameter X and the jacobian determinant in dX
  void nodalShape(double X[], double &dX, double xi, double eta);
  void nodalShape(arma::vec &X, double &dX, double xi, double eta);
  void nodalShape(arma::vec &X, arma::mat &jac,  double &dX, double xi, double eta);

  // compute the b-moments using the values already assigned in the object
  void compute_moments();

  // compute the b-moments for the specified f function
  void compute_moments(double (*f)(double, double));

  // compute the b-moments for the Fval function values
  void compute_moments(arma::vec Fval);
};

// Future class definition
class BMoment3D
{
};

#endif