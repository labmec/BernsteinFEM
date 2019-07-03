#include "Moments.h"


// constructors
template <typename _s, Element_t EL>
BMoment<_s, EL>::BMoment(uint q, uint n, const Element<EL> &element)
    : Bmoment(), Cval(), element(element), f()
{
    this->q = q;
    this->n = n;
    this->nb_Array = nb_Array;
    this->element.setPermutationPOrder(n);
}

// copy constructor
template <typename _s, Element_t EL>
BMoment<_s, EL>::BMoment(const BMoment<_s, EL> &cp)
    : Bmoment(cp.Bmoment), Cval(cp.Cval), element(cp.element), f(cp.f)
{
    this->q = cp.q;
    this->n = cp.n;
    this->nb_Array = cp.nb_Array;
    this->element.setPermutationPOrder(n);
}

// copy asignment operator
template <typename _s, Element_t EL>
BMoment<_s, EL> &BMoment<_s, EL>::operator=(const BMoment &cp)
{
    if (this != &cp)
    {
        this->q = cp.q;
        this->n = cp.n;
        this->nb_Array = cp.nb_Array;
    }
    return *this;
}

// destructor
template <typename _s, Element_t EL>
BMoment<_s, EL>::~BMoment() {}

// getters
// returns number of integration points
template <typename _s, Element_t EL>
uint BMoment<_s, EL>::getNumIntegrationPoints() { return q; }

// returns polynomial order
template <typename _s, Element_t EL>
uint BMoment<_s, EL>::getPOrder() { return n; }

// returns moments array length
template <typename _s, Element_t EL>
uint BMoment<_s, EL>::getLenMoments() { return lenMoments; }

// returns number
template <typename _s, Element_t EL>
uint BMoment<_s, EL>::getLenCval() { return lenCval; }

// return the dimension of function Image (function is scalar valued by default)
template <typename _s, Element_t EL>
uint BMoment<_s, EL>::getNbArray() { return nb_Array; }

// returns wether you will be using function values or function definition (true for function values)
template <typename _s, Element_t EL>
bool BMoment<_s, EL>::getFunctVal() { return functVal; }

// returns the element used for computing
template <typename _s, Element_t EL>
Element<EL> &BMoment<_s, EL>::getElement() { return element; }

// returns a copy of the Bmoment matrix
template <typename _s, Element_t EL>
const TPZVec<REAL> &BMoment<_s, EL>::getMoments() { return Bmoment; }

// setters
// sets number of integration points
template <typename _s, Element_t EL>
void BMoment<_s, EL>::setNumIntegrationPoints(uint q) { this->q = q; }

// sets polynomial order
template <typename _s, Element_t EL>
void BMoment<_s, EL>::setPOder(uint n) { this->n = n; }

// sets the dimension of function Image (function is scalar valued by default)
template <typename _s, Element_t EL>
void BMoment<_s, EL>::setNbArray(uint nb_Array)
{
    this->nb_Array = nb_Array;
    Bmoment.resize(lenMoments, nb_Array);
    Cval.resize(lenCval, nb_Array);
}

template <typename _s, Element_t EL>
void BMoment<_s, EL>::setFunctionValues(const TPZVec<REAL> &Cval)
{
    this->Cval = Cval;
    fValSet = true;
}

template <typename _s, Element_t EL>
void BMoment<_s, EL>::setFunctionDefinition(std::function<_s> f) { this->f = f; }

template <typename _s, Element_t EL>
void BMoment<_s, EL>::setElement(Element<EL> element) { this->element = element; }

template <typename _s, Element_t EL>
void BMoment<_s, EL>::zero() { Bmoment.zeros(); }

template <typename _s, Element_t EL>
void BMoment<_s, EL>::useFunctionDef() { functVal = false; }

template <typename _s, Element_t EL>
void BMoment<_s, EL>::useFunctionValue() { functVal = true; }

template <typename _s, Element_t EL>
void BMoment<_s, EL>::computeMoments(std::function<_s> f)
{
    setFunctionDefinition(f);
    computeMoments();
}

template <typename _s, Element_t EL>
void BMoment<_s, EL>::computeMoments(const TPZVec<REAL> &Cval)
{
    setFunctionValues(Cval);
    computeMoments();
}

// instantiate the types that will be used throughout the library
template class BMoment<double(double), Element_t::LinearEl>;
template class BMoment<double (double, double), Element_t::TriangularEl>;
template class BMoment<double (double, double), Element_t::QuadrilateralEl>;
template class BMoment<double(double, double, double), Element_t::CubeEl>;
template class BMoment<double(double, double, double), Element_t::TetrahedronEl>;