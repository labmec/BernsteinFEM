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
#include "Common/pzreal.h"
#include <vector>
#include "Element_t.h"

#ifndef uint
using uint = unsigned;
#endif // !uint

template <Element_t El>
class Permutation
{
private:
  uint n;                           // polynomial order, in order to compute the vectors
  std::vector<uint> permutationVec; // permutation vector to organize elements of the matrices
  std::vector<uint> inversePermVec; // inverse permutation vector
  std::vector<uint> *idxVec;        // index vector for each vertex (the same as in the element which holds this object)
  bool idxVecInstantiated;          // indicates whether we instantiated a new idxVec or just referenced an object
  bool pVecComputed = false;        // indicates whether the permutation vector is computed or not
  bool i_pVecComputed = false;      // indicates whether the inverse permutation vector is computed or not

public:
  Permutation();

  Permutation(const uint &pOrder, std::vector<uint> &idxVec);

  Permutation(std::vector<uint> &idxVec);

  Permutation(Permutation const &cp);

  ~Permutation();

  // getters

  uint getPOrder();
  std::vector<uint> &getIndexVector();
  std::vector<uint> &getPermutationVector();
  std::vector<uint> &getInvPermutationVec();

  //setters
  void setPOrder(uint n);

  // sets the index vector of the object, according to the indices of the vertices of an element
  void setIndexVector(std::vector<uint> &idxVector);

  // computes the permutation vector
  void computePermVec();

  // computes the inverse permutation vector
  void computeInvPermVec();
};