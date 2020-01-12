/***************************************************************
 * Author: Lucas B. Andrade
 * 
 * This file contains definitions of a Elements template class
 * The template argument is an enum Element_t, which specifies
 * the type of element of the object
 * 
 * The types of elements are in respect to its geometrical
 * properties (e.g. LinearEl is the 1D line element)
 ***************************************************************/

#pragma once

#include "pzreal.h"     // REAL, STATE
#include "pzmatrix.h"   // TPZFMatrix
#include "pzvec.h"      // TPZVec
#include "Permutations.h"
#include "Element_t.h"

#ifndef uint
using uint = unsigned;
#endif

// simplicial elements (e.g. triangle, tetrahedron) are mapped
// using the Duffy Transform, from the equivalent quadrilateral/cube master element
// e.g. you should use the square [0,1] X [0,1] as master element to the TriangleEl

template <Element_t EL>
class Element
{
    TPZFMatrix<REAL> vertices;    // vertices of the element
    TPZFMatrix<REAL> coordinates; // coordinates object to return from 'mapElement'
    TPZVec<int64_t> idxVec;      // stores index for each vertex
    static TPZFMatrix<REAL> jac;  // matrix to store the values of the jacobian when no arguments are passed
    Permutation<EL> perm;         // permutation vector for positioning

public:
    // constructors
    // default constructor, makes a element equivalent to the Master Element
    Element();

    // constructor by the vertices of the element
    Element(TPZFMatrix<REAL> &vertices);

    // constructor by the vertices and indices of the vertices of the element
    Element(TPZFMatrix<REAL> &vertices, TPZVec<int64_t> &indexVector) : vertices(vertices), idxVec(indexVector), perm(indexVector)
    {
    }

    // copy constructor
    Element(const Element<EL> &cp);

    // copy assignment operator
    Element<EL> &operator=(const Element<EL> &cp)
    {
        if (this != &cp)
            vertices = cp.vertices;
        return *this;
    }

	~Element() = default;

    TPZVec<int64_t> &getIndexVector() { return idxVec; }

    void setPermutationPOrder(uint64_t n)
    {
        perm.setPOrder(n);
        perm.computePermVec();
    }

    void setIndexVector(TPZVec<int64_t> &idxVec)
    {
        this->idxVec = idxVec;
        perm.setIndexVector(idxVec);
    }

    // returns the vertices of the element object
    TPZFMatrix<REAL> &getVertices() { return vertices; }

    // get last jacobian matrix that was computed
    static TPZFMatrix<REAL> &getLastJacobian() { return jac; }

    // maps the specified point of the element
    // to the position in which it should be
    // in the vectors and matrices throughout the computations
    uint64_t position(const TPZVec<uint64_t> &point);

	// same as above but with initializer_list
	uint64_t position(const std::initializer_list<uint64_t> &point);

    // Maps the element from the matrix element
    // Given the coordinates in the master element in [0,1],
    // returns the corresponding coordinates in the element
    TPZFMatrix<REAL> &mapToElement(TPZFMatrix<REAL> &coordinates, TPZFMatrix<REAL> &jacobian = jac);

	// same as above but with initializer list
	TPZFMatrix<REAL> &mapToElement(const std::initializer_list<REAL> &coordinates, TPZFMatrix<REAL> &jacobian = jac);
};