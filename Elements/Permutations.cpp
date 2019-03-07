#include "Permutations.h"

// generic implementations

template <Element_t El> 
Permutation<El>::Permutation() { }

template <Element_t El>
Permutation<El>::Permutation(arma::ivec &idxVec) : idxVec(&idxVec) { }

template <Element_t El>
Permutation<El>::Permutation(Permutation const &cp)
{
    this->n = cp.n;
}

template <Element_t El>
uint Permutation<El>::getPOrder()
{
    return n;
}

template <Element_t El>
std::vector<uint> &Permutation<El>::getPermutationVector()
{
    if (!pVecComputed)
        computePermVec();

    return permutationVec;
}

template <Element_t El>
std::vector<uint> &Permutation<El>::getInvPermutationVec()
{
    if (!i_pVecComputed)
        computeInvPermVec();

    return inversePermVec;
}

template <Element_t El>
void Permutation<El>::setPOrder(uint n)
{
    if (this->n != n)
    {
        this->n = n;
        pVecComputed = false;
        i_pVecComputed = false;
    }
}

template <Element_t El>
void Permutation<El>::setIndexVector(arma::ivec &idxVec)
{
    this->idxVec = &idxVec;
}

// template class instantiations

template class Permutation<Element_t::LinearEl>;
template class Permutation<Element_t::TriangularEl>;
template class Permutation<Element_t::QuadrilateralEl>;
template class Permutation<Element_t::TetrahedronEl>;
template class Permutation<Element_t::CubeEl>;

// specializations

// LinearEl

template <>
void Permutation<Element_t::LinearEl>::computePermVec()
{
    permutationVec.resize(n);

    std::iota(permutationVec.begin(), permutationVec.end(), 0);
}

template <>
void Permutation<Element_t::LinearEl>::computeInvPermVec()
{
    inversePermVec.resize(n);

    std::iota(inversePermVec.begin(), inversePermVec.end(), 0);
}

// TriangularEl

template <>
void Permutation<Element_t::TriangularEl>::computePermVec()
{
    permutationVec.resize((n + 1) * (n + 1));

    uint elAdded = 0;

    permutationVec[0] = 0;           // v0
    permutationVec[n * (n + 1)] = 1; // v1
    permutationVec[n] = 2;           // v2

    elAdded += 3;

    // L0
    // TODO: check orientation
    if (idxVec->at(0) < idxVec->at(1))
    {
        for (uint a1 = 1; a1 < n; a1++)
        {
            permutationVec[a1 * (n + 1)] = elAdded++;
        }
    } else
    {
        for (uint a1 = n - 1; a1 > 0; a1--)
        {
            permutationVec[a1 * (n + 1)] = elAdded++;
        }
    }
    

    // L1
    // TODO: check orientation
    if (idxVec->at(1) < idxVec->at(2))
    {
        for (uint a1 = 1; a1 < n; a1++)
        {
            uint a2 = n - a1;
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    } else
    {
        for (uint a1 = n - 1; a1 > 0; a1--)
        {
            uint a2 = n - a1;
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    }
    

    // L2
    // TODO: check orientation
    if (idxVec->at(2) < idxVec->at(0))
    {
        for (uint a2 = 1; a2 < n; a2++)
        {
            permutationVec[a2] = elAdded++;
        }
    } else
    {
        for (uint a2 = n - 1; a2 > 0; a2--)
        {
            permutationVec[a2] = elAdded++;
        }
    }
    

    // Middle points
    for (uint a1 = 1; a1 < n; a1++)
    {
        for (uint a2 = 1; a2 < n - a1; a2++)
        {
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    }

    // std::cout << "Number of elements that were not computed for permutation: " << ((n + 1) * (n + 2) * 0.5) - elAdded << std::endl;

    // std::cout << "Permutation vector:" << std::endl;
    // for (uint i : permutationVec)
    // {
    //     std::cout << i << " ";
    // }
    // std::cout << std::endl;
}

template <>
void Permutation<Element_t::TriangularEl>::computeInvPermVec()
{
    inversePermVec.resize((n + 1) * (n + 2) * 0.5);

    // uint elAdded = 0;
}

// QuadrilateralEl

template <>
void Permutation<Element_t::QuadrilateralEl>::computePermVec()
{
    permutationVec.resize((n + 1) * (n + 1));

    int elAdded = 0;

    permutationVec[0] = 0;               // v0
    permutationVec[n] = 1;               // v1
    permutationVec[(n + 1) * n + n] = 2; // v2
    permutationVec[(n + 1) * n] = 3;     // v3

    elAdded += 4;

    // L0
    // TODO: check orientation
    if (idxVec->at(0) < idxVec->at(1))
    {
        for (uint i = 1; i < n; i++)
        {
            // permutationVec[i + 1] = 4 + i;
            permutationVec[i] = elAdded++;
        }
    } else
    {
        for (uint i = n - 1; i > 0; i--)
        {
            permutationVec[i] = elAdded++;
        }
    }
    

    // L1
    // TODO: check orientation
    if (idxVec->at(1) < idxVec->at(2))
    {
        for (uint i = 2; i <= n; i++)
        {
            permutationVec[i * (n + 1) - 1] = elAdded++;
        }
    } else
    {
        for (uint i = n; i > 2; i--)
        {
            permutationVec[i * (n + 1) - 1] = elAdded++;
        }
    }
    

    // L2
    // TODO: check orientation
    uint ini = (n + 1) * n;
    if (idxVec->at(2) >= idxVec->at(3))
    {
        for (uint i = 1; i < n; i++)
        {
            permutationVec[ini + i] = elAdded++;
        }
    } else
    {
        for (uint i = n - 1; i > 0; i--)
        {
            permutationVec[ini + i] = elAdded++;
        }
    }
    

    // L3
    // TODO: check orientation
    if (idxVec->at(3) >= idxVec->at(0))
    {
        for (uint i = 1; i < n; i++)
        {
            permutationVec[i * (n + 1)] = elAdded++;
        }
    } else
    {
        for (uint i = n - 1; i > 0; i--)
        {
            permutationVec[i * (n + 1)] = elAdded++;
        }
    }
    

    // Middle points
    for (uint a1 = 1; a1 < n; a1++)
    {
        for (uint a2 = 1; a2 < n; a2++)
        {
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    }

    // std::cout << "Number of elements that were not computed for permutation: " << permutationVec.size() - elAdded << std::endl;

    // std::cout << "Permutation vector:" << std::endl;
    // for (uint i : permutationVec)
    // {
    //     std::cout << i << " ";
    // }
    // std::cout << std::endl;
}

template <>
void Permutation<Element_t::QuadrilateralEl>::computeInvPermVec()
{
    inversePermVec.resize((n + 1) * (n + 1));

    inversePermVec[0] = 0;
    inversePermVec[1] = n;
    inversePermVec[2] = (n + 1) * n + n;
    inversePermVec[3] = (n + 1) * n;
    // NOTE: not implemented because we are doing permutation in-place
}

template <>
void Permutation<Element_t::CubeEl>::computePermVec()
{

}

template <>
void Permutation<Element_t::CubeEl>::computeInvPermVec()
{
}

template <>
void Permutation<Element_t::TetrahedronEl>::computePermVec()
{
}

template <>
void Permutation<Element_t::TetrahedronEl>::computeInvPermVec()
{
}