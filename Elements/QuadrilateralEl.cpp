#include "Elements.h"

#define QEL Element_t::QuadrilateralEl

template <>
TPZFMatrix<REAL> Element<QEL>::jac(2, 2, 0);

template <>
Element<QEL>::Element()
    : vertices(4, 2),
    coordinates(1, 2),
	idxVec({ 0, 1, 2, 3 }),
    perm(idxVec)
{
    vertices(0, 0) = 0.0; vertices(0, 1) = 0.0; // v1 = (0,0)
    vertices(1, 0) = 1.0; vertices(1, 1) = 0.0; // v2 = (1,0)
    vertices(2, 0) = 1.0; vertices(2, 1) = 1.0; // v3 = (1,1)
    vertices(3, 0) = 0.0; vertices(3, 1) = 1.0; // v4 = (0,1)
}

template <>
Element<QEL>::Element(TPZFMatrix<REAL> &v)
    : vertices(4, 2),
    coordinates(1, 2),
	idxVec({ 0, 1, 2, 3 }),
    perm(idxVec)
{
    if (v.Rows() < 4 || v.Cols() < 2)
    {
        throw std::invalid_argument("QuadrilateralEl vertices constructor: not enough size in argument (at least 4x2)");
    }
    else
    {
        vertices(0, 0) = v(0, 0);
        vertices(0, 1) = v(0, 1);
        vertices(1, 0) = v(1, 0);
        vertices(1, 1) = v(1, 1);
        vertices(2, 0) = v(2, 0);
        vertices(2, 1) = v(2, 1);
        vertices(3, 0) = v(3, 0);
        vertices(3, 1) = v(3, 1);
    }
}

template <>
Element<QEL>::Element(const Element<QEL> &cp)
	: vertices(cp.vertices),
	coordinates(1, 2),
	idxVec(cp.idxVec),
	perm(cp.perm)
{}

template <>
uint64_t Element<QEL>::position(const TPZVec<uint64_t> &point)
{
    uint64_t n = perm.getPOrder();
    if (point.size() >= 2)
        return perm.getPermutationVector()[ point[0] * (n + 1) + point[1] ];
    else
        throw std::logic_error("Quadrilateral Element 'position' method called with too few vector elements\n\t2 required");
}

uint64_t Element<QEL>::position(const std::initializer_list<uint64_t> &point) {
	uint64_t n = perm.getPOrder();
	if (point.size() >= 2) {
		auto it = point.begin();
		const uint64_t a1 = *it, a2 = *(++it);
		return perm.getInvPermutationVec()[ a1 * (n + 1) + a2 ];
	} else {
		throw std::logic_error("Quadrilateral Element 'position' method called with too few vector elements\n\t2 required");
	}
}

// implemented using nodal shape function
template <>
TPZFMatrix<REAL> &Element<QEL>::mapToElement(TPZFMatrix<REAL> &xi, TPZFMatrix<REAL> &jacobian)
{
    double xi_ = xi(0);
    double eta_ = xi(1);
    {
        double N0 = (1.0 - xi_) * (1.0 - eta_);
        double N1 = xi_ * (1.0 - eta_);
        double N2 = xi_ * eta_;
        double N3 = (1.0 - xi_) * eta_;

        coordinates(0) = N0 * vertices(0, 0) + N1 * vertices(1, 0) + N2 * vertices(2, 0) + N3 * vertices(3, 0);
        coordinates(1) = N0 * vertices(0, 1) + N1 * vertices(1, 1) + N2 * vertices(2, 1) + N3 * vertices(3, 1);
    }

    jacobian(0, 0) = (1.0 - eta_) * (vertices(1, 0) - vertices(0, 0)) + eta_ * (vertices(2, 0) - vertices(3, 0));
    jacobian(0, 1) = (1.0 - xi_) * (vertices(3, 0) - vertices(0, 0)) + xi_ * (vertices(2, 0) - vertices(1, 0));
    jacobian(1, 0) = (1.0 - eta_) * (vertices(1, 1) - vertices(0, 1)) + eta_ * (vertices(2, 1) - vertices(3, 1));
    jacobian(1, 1) = (1.0 - xi_) * (vertices(3, 1) - vertices(0, 1)) + xi_ * (vertices(2, 1) - vertices(1, 1));

    return coordinates;
}

template <>
TPZFMatrix<REAL>& Element<QEL>::mapToElement(const std::initializer_list<REAL>& xi, TPZFMatrix<REAL>& jacobian)
{
	auto it = xi.begin();

	double xi_ = *it;
	double eta_ = *(++it);
	{
		double N0 = (1.0 - xi_) * (1.0 - eta_);
		double N1 = xi_ * (1.0 - eta_);
		double N2 = xi_ * eta_;
		double N3 = (1.0 - xi_) * eta_;

		coordinates(0) = N0 * vertices(0, 0) + N1 * vertices(1, 0) + N2 * vertices(2, 0) + N3 * vertices(3, 0);
		coordinates(1) = N0 * vertices(0, 1) + N1 * vertices(1, 1) + N2 * vertices(2, 1) + N3 * vertices(3, 1);
	}

	jacobian(0, 0) = (1.0 - eta_) * (vertices(1, 0) - vertices(0, 0)) + eta_ * (vertices(2, 0) - vertices(3, 0));
	jacobian(0, 1) = (1.0 - xi_) * (vertices(3, 0) - vertices(0, 0)) + xi_ * (vertices(2, 0) - vertices(1, 0));
	jacobian(1, 0) = (1.0 - eta_) * (vertices(1, 1) - vertices(0, 1)) + eta_ * (vertices(2, 1) - vertices(3, 1));
	jacobian(1, 1) = (1.0 - xi_) * (vertices(3, 1) - vertices(0, 1)) + xi_ * (vertices(2, 1) - vertices(1, 1));

	return coordinates;
}