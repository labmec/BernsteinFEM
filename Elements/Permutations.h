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

#include <armadillo>
#include <vector>
#include "Element_t.h"

#ifndef uint
#define uint unsigned
#endif

template <Element_t El>
class Permutation
{
private:
  uint n;                           // polynomial order, in order to compute the vectors
  std::vector<uint> permutationVec; // permutation vector to organize elements of the matrices
  std::vector<uint> inversePermVec; // inverse permutation vector
  arma::ivec *idxVec;                // index vector for each vertex (the same as in the element which holds this object)
  bool pVecComputed = false;        // indicates whether the permutation vector is computed or not
  bool i_pVecComputed = false;      // indicates whether the inverse permutation vector is computed or not

  void computePermVec();
  void computeInvPermVec();

public:
  Permutation();

  Permutation(const uint &n);

  Permutation(arma::ivec &idxVec);

  Permutation(Permutation const &cp);

  // getters

  uint getPOrder();
  std::vector<uint> &getPermutationVector();
  std::vector<uint> &getInvPermutationVec();

  //setters
  void setPOrder(uint n);

  // sets the index vector of the object, according to the indices of the vertices of an element
  void setIndexVector(arma::ivec &idxVector);
};

template <Element_t El>
class PermutationPool
{
  private:
    static std::vector<Permutation<El> *> pool;

  public:
    static Permutation<El> &GetPermutation(uint n)
    {
      if (n <= 1)
      {
        throw std::logic_error("GetPermutation(n) only accepts n > 1");
      }
      if (pool.size() > n - 1)
      {
        if (pool[n - 2] == nullptr)
        {
          return *(pool[n - 2] = new Permutation<El>(n));
        } else {
          return *(pool[n - 2]);
        }
      } else {
        pool.resize(pool.size() * 2);
        return *(pool[n - 2] = new Permutation<El>(n));
      }
    }
};