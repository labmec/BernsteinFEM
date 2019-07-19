/***************************************************************
 * Author: Lucas B. Andrade
 * 
 * This file contains the definition of the Permutation class
 * which holds permutation vectors using indices for permutation
 * complying with organization of the element matrices used throughout
 * the NeoPZ environment
 * 
 * The permutation class is used exclusively as an object in the Element class
 * making the permutation in-place through the 'position' method
 ***************************************************************/

#pragma once

#define REALdouble
#include "pzreal.h"
#include "pzvec.h"
#include "Element_t.h"

template <Element_t El>
class Permutation
{
private:
  uint64_t n;                       // polynomial order, in order to compute the vectors
  TPZVec<uint64_t> permutationVec;  // permutation vector to organize elements of the matrices
  TPZVec<uint64_t> inversePermVec;  // inverse permutation vector
  TPZVec<uint64_t> *idxVec;         // index vector for each vertex (the same as in the element which holds this object)
  bool idxVecInstantiated;          // indicates whether we instantiated a new idxVec or just referenced an object
  bool pVecComputed = false;        // indicates whether the permutation vector is computed or not
  bool i_pVecComputed = false;      // indicates whether the inverse permutation vector is computed or not

public:
  Permutation();

  Permutation(const uint64_t &pOrder, TPZVec<uint64_t> &idxVec);

  Permutation(TPZVec<uint64_t> &idxVec);

  Permutation(Permutation const &cp);

  ~Permutation();

  // getters

  uint64_t getPOrder();
  TPZVec<uint64_t> &getIndexVector();
  TPZVec<uint64_t> &getPermutationVector();
  TPZVec<uint64_t> &getInvPermutationVec();

  //setters
  void setPOrder(uint64_t n);

  // sets the index vector of the object, according to the indices of the vertices of an element
  void setIndexVector(TPZVec<uint64_t> &idxVector);

  // computes the permutation vector
  void computePermVec();

  // computes the inverse permutation vector
  void computeInvPermVec();
};