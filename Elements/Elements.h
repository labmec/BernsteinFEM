#pragma once
#include <armadillo>

enum Element_t
{
  LinearEl,
  QuadrilateralEl,
  TriangularEl,
  CubeEl,
  TetrahedronEl
};

// simplicial elements (e.g. triangle, tetrahedron) are mapped
// using the Duffy Transform, from the equivalent quadrilateral/cube master element
// e.g. you should use the square [0,1] X [0,1] as master element to the TriangleEl

template <Element_t EL>
class Element
{
  arma::mat vertices;    // vertices of the element
  arma::mat coordinates; // coordinates object to return from 'mapElement'
  static arma::mat jac;  // matrix to store the values of the jacobian when no arguments are passed

public:
  // constructors
  // default constructor, makes a element equal to the Master Element
  Element();

  // constructor by the vertices of the element, the first coord
  Element(const arma::mat &vertices);

  // copy constructor
  Element(const Element<EL> &cp);

  // copy assignment operator
  Element &operator=(const Element& cp)
  {
    if (this != &cp)
    vertices = cp.vertices;
    return *this;
  }

  // get last jacobian matrix that was computed
  static const arma::mat &getLastJacobian() { return jac; }

  // Maps the element from the matrix element
  // Given the coordinates in the master element in [0,1],
  // returns the corresponding coordinates in the element
  const arma::mat &mapToElement(const arma::mat &coordinates, arma::mat &jacobian = jac);
};