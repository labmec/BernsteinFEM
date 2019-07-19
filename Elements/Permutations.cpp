#include "Permutations.h"

// generic implementations

template <Element_t El> 
Permutation<El>::Permutation() 
{
	this->idxVec = new TPZVec<uint64_t>({ 0,1,2,3 }); 
	idxVecInstantiated = true;
}

template <Element_t El>
Permutation<El>::Permutation(const uint64_t &n, TPZVec<uint64_t> &idxVec) : n(n), idxVec(&idxVec)
{
    computePermVec();
	idxVecInstantiated = false;
}

template <Element_t El>
Permutation<El>::Permutation(TPZVec<uint64_t> &idxVec) : idxVec(&idxVec), idxVecInstantiated(false) { }

template <Element_t El>
Permutation<El>::Permutation(Permutation const &cp) : n(cp.n), idxVec(cp.idxVec), idxVecInstantiated(false) { }

template <Element_t El>
Permutation<El>::~Permutation()
{
	if (idxVecInstantiated)
	{
		free(idxVec);
	}
}

template <Element_t El>
uint64_t Permutation<El>::getPOrder()
{
    return n;
}

template <Element_t El>
TPZVec<uint64_t> &Permutation<El>::getIndexVector()
{
    return *idxVec;
}

template <Element_t El>
TPZVec<uint64_t> &Permutation<El>::getPermutationVector()
{
    return permutationVec;
}

template <Element_t El>
TPZVec<uint64_t> &Permutation<El>::getInvPermutationVec()
{
    if (!i_pVecComputed)
        computeInvPermVec();

    return inversePermVec;
}

template <Element_t El>
void Permutation<El>::setPOrder(uint64_t n)
{
    if (this->n != n)
    {
        this->n = n;
		computePermVec();
        i_pVecComputed = false;
    }
}

template <Element_t El>
void Permutation<El>::setIndexVector(TPZVec<uint64_t> &idxVec)
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
    permutationVec.resize(n + 1);

	for (unsigned i = 0; i < permutationVec.size(); permutationVec[i] = i++);

	pVecComputed = true;
}

template <>
void Permutation<Element_t::LinearEl>::computeInvPermVec()
{
    inversePermVec.resize(n + 1);

	for (unsigned i = 0; i < permutationVec.size(); permutationVec[i] = i++);

	i_pVecComputed = true;
}

// TriangularEl

template <>
void Permutation<Element_t::TriangularEl>::computePermVec()
{
    permutationVec.resize((n + 1) * (n + 1));

    uint64_t elAdded = 0;

    permutationVec[0] = 0;           // v0
    permutationVec[n * (n + 1)] = 1; // v1
    permutationVec[n] = 2;           // v2

    elAdded += 3;

    // L0
    // TODO: check orientation
    if (idxVec->operator[](0) < idxVec->operator[](1))
    {
        for (uint64_t a1 = 1; a1 < n; a1++)
        {
            permutationVec[a1 * (n + 1)] = elAdded++;
        }
    } else {
        for (uint64_t a1 = n - 1; a1 > 0; a1--)
        {
            permutationVec[a1 * (n + 1)] = elAdded++;
        }
    }
    

    // L1
    // TODO: check orientation
    if (idxVec->operator[](1) < idxVec->operator[](2))
    {
        for (uint64_t a1 = 1; a1 < n; a1++)
        {
            uint64_t a2 = n - a1;
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    } else {
        for (uint64_t a1 = n - 1; a1 > 0; a1--)
        {
            uint64_t a2 = n - a1;
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    }
    

    // L2
    // TODO: check orientation
    if (idxVec->operator[](2) < idxVec->operator[](0))
    {
        for (uint64_t a2 = 1; a2 < n; a2++)
        {
            permutationVec[a2] = elAdded++;
        }
    } else {
        for (uint64_t a2 = n - 1; a2 > 0; a2--)
        {
            permutationVec[a2] = elAdded++;
        }
    }
    

    // Middle points
    for (uint64_t a1 = 1; a1 < n; a1++)
    {
        for (uint64_t a2 = 1; a2 < n - a1; a2++)
        {
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    }

	pVecComputed = true;

    // std::cout << "Number of elements that were not computed for permutation: " << ((n + 1) * (n + 2) * 0.5) - elAdded << std::endl;

    // std::cout << "Permutation vector:" << std::endl;
    // for (uint64_t i : permutationVec)
    // {
    //     std::cout << i << " ";
    // }
    // std::cout << std::endl;
}

template <>
void Permutation<Element_t::TriangularEl>::computeInvPermVec()
{
    inversePermVec.resize((n + 1) * (n + 2) * 0.5);

    // uint64_t elAdded = 0;
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
    if (idxVec->operator[](0) < idxVec->operator[](1))
    {
        for (uint64_t i = 1; i < n; i++)
        {
            // permutationVec[i + 1] = 4 + i;
            permutationVec[i] = elAdded++;
        }
    } else {
        for (uint64_t i = n - 1; i > 0; i--)
        {
            permutationVec[i] = elAdded++;
        }
    }
    

    // L1
    // TODO: check orientation
    if (idxVec->operator[](1) < idxVec->operator[](2))
    {
        for (uint64_t i = 2; i <= n; i++)
        {
            permutationVec[i * (n + 1) - 1] = elAdded++;
        }
    } else {
        for (uint64_t i = n; i > 2; i--)
        {
            permutationVec[i * (n + 1) - 1] = elAdded++;
        }
    }
    

    // L2
    // TODO: check orientation
    uint64_t ini = (n + 1) * n;
    if (idxVec->operator[](2) >= idxVec->operator[](3))
    {
        for (uint64_t i = 1; i < n; i++)
        {
            permutationVec[ini + i] = elAdded++;
        }
    } else {
        for (uint64_t i = n - 1; i > 0; i--)
        {
            permutationVec[ini + i] = elAdded++;
        }
    }
    

    // L3
    // TODO: check orientation
    if (idxVec->operator[](3) >= idxVec->operator[](0))
    {
        for (uint64_t i = 1; i < n; i++)
        {
            permutationVec[i * (n + 1)] = elAdded++;
        }
    } else {
        for (uint64_t i = n - 1; i > 0; i--)
        {
            permutationVec[i * (n + 1)] = elAdded++;
        }
    }
    

    // Middle points
    for (uint64_t a1 = 1; a1 < n; a1++)
    {
        for (uint64_t a2 = 1; a2 < n; a2++)
        {
            permutationVec[a1 * (n + 1) + a2] = elAdded++;
        }
    }

	pVecComputed = true;

    // std::cout << "Number of elements that were not computed for permutation: " << permutationVec.size() - elAdded << std::endl;

    // std::cout << "Permutation vector:" << std::endl;
    // for (uint64_t i : permutationVec)
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
	// NOTE2: might be faster to not do in-place, although
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