#include "Moments.h"

// constructors
template <typename _s, Element_t EL>
BMoment<_s, EL>::BMoment(int q, int n, Element<EL> element = Element<EL>(), int nb_Array = 1)
    : Bmoment(), Cval(), element(element), quadraWN(), f()
{
    this->q = q;
    this->n = n;
    this->nb_Array = nb_Array;
}

// copy constructor
template <typename _s, Element_t EL>
BMoment<_s, EL>::BMoment(const BMoment<double(double), Element_t::LinearEl> &cp)
{
    this->q = cp.q;
    this->n = cp.n;
    this->nb_Array = cp.nb_Array;
}

// destructor
template <typename _s, Element_t EL>
BMoment<_s, EL>::~BMoment() {}

// getters
// returns number of integration points
template <typename _s, Element_t EL>
inline int BMoment<_s, EL>::getNumIntegrationPoints() { return q; }

// returns polynomial order
template <typename _s, Element_t EL>
inline int BMoment<_s, EL>::getPOrder() { return n; }

// returns moments array length
template <typename _s, Element_t EL>
inline int BMoment<_s, EL>::getLenMoments() { return lenMoments; }

// returns number
template <typename _s, Element_t EL>
inline int BMoment<_s, EL>::getLenCval() { return lenCval; }

// return the dimension of function Image (function is scalar valued by default)
template <typename _s, Element_t EL>
inline int BMoment<_s, EL>::getNbArray() { return nb_Array; }

// returns wether you will be using function values or function definition (true for function values)
template <typename _s, Element_t EL>
inline bool BMoment<_s, EL>::getFunctVal() { return functVal; }

// returns the element used for computing
template <typename _s, Element_t EL>
inline Element<EL> BMoment<_s, EL>::getElement() { return element; }

// returns the whole Bmoment matrix
template <typename _s, Element_t EL>
inline const arma::mat &BMoment<_s, EL>::getBMoment() { return Bmoment; }

// setters
// sets number of integration points
template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::setNumIntegrationPoints(int q) { this->q = q; }

// sets polynomial order
template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::setPOder(int n) { this->n = n; }

// sets the dimension of function Image (function is scalar valued by default)
template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::setNbArray(int nb_Array)
{
    this->nb_Array = nb_Array;
    Bmoment.resize(lenMoments, nb_Array);
    Cval.resize(lenCval, nb_Array);
}

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::setFunctionValues(const arma::vec &Cval) { this->Cval.swap(Cval); }

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::setFunctionValues(const arma::mat &Cval) { this->Cval.swap(Cval); }

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::setFunctionDefinition(std::function<_s> f) { this->f = f; }

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::setElement(Element<EL> element) { this->element = element; }

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::zero() { Bmoment.zeros(); }

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::useFunctionDef() { functVal = false; }

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::useFunctionValue() { functVal = true; }

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::computeMoments(std::function<_s> f)
{
    setFunctionDefinition(f);
    compute_moments();
}

template <typename _s, Element_t EL>
inline void BMoment<_s, EL>::computeMoments(const arma::mat &Cval)
{
    setFunctionValues(Cval);
    compute_moments();
}

template <typename _s, Element_t EL>
inline int position(int i)
{
    return position(i, this->n);
}