#include "Evaluation.h"

template <Element_t EL>
BEval<EL>::BEval(uint q, uint n, const Element<EL> &el)
    : BBVec(), eval(), element(el) { }

template <Element_t EL>
BEval<EL>::BEval(uint q, uint n, const arma::vec &coeffVec, const Element<EL> &el)
    : BBVec(coeffVec), eval(), element(el) { }

template <Element_t EL>
uint BEval<EL>::getNumIntegrationPoints()
{
    return q;
}

template <Element_t EL>
uint BEval<EL>::getPOrder()
{
    return n;
}

template <Element_t EL>
arma::vec &BEval<EL>::getCoefficientsVec()
{
    return BBVec;
}

template <Element_t EL>
arma::vec &BEval<EL>::getEval()
{
    return eval;
}

template <Element_t EL>
Element<EL> &BEval<EL>::getElement()
{
    return element;
}

template <Element_t EL>
void BEval<EL>::setCoefficientsVec(const arma::vec &vec)
{
    if (vec.size() >= bbvec_len)
        BBVec = vec;
    else
    {
        std::stringstream str;
        str << "Coefficient vector passed for 'BEval' object should have " 
            << bbvec_len << " elements, has: " << vec.size() << "\n";
        throw new std::logic_error(str.str().c_str());
    }
}

// instantiate the types that will be used throughout the library
template class BEval<Element_t::LinearEl>;
template class BEval<Element_t::TriangularEl>;
template class BEval<Element_t::QuadrilateralEl>;
template class BEval<Element_t::CubeEl>;
template class BEval<Element_t::TetrahedronEl>;