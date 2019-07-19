#include "Elements.h"

#include "Elements.h"

#define CEL Element_t::CubeEl

template <>
TPZFMatrix<REAL> Element<CEL>::jac(3, 3, 0);

template <>
Element<CEL>::Element()
	: vertices(8, 3),
	coordinates(1, 3),
	idxVec({ 0,1,2,3,4,5 }),
	perm(idxVec)
{
    vertices(0, 0) = 0.0; vertices(0, 1) = 0.0; vertices(0, 2) = 0.0;
    vertices(1, 0) = 1.0; vertices(1, 1) = 0.0; vertices(1, 2) = 0.0;
    vertices(2, 0) = 1.0; vertices(2, 1) = 1.0; vertices(2, 2) = 0.0;
    vertices(3, 0) = 0.0; vertices(3, 1) = 1.0; vertices(3, 2) = 0.0;
    vertices(4, 0) = 0.0; vertices(4, 1) = 0.0; vertices(4, 2) = 1.0;
    vertices(5, 0) = 1.0; vertices(5, 1) = 0.0; vertices(5, 2) = 1.0;
    vertices(6, 0) = 1.0; vertices(6, 1) = 1.0; vertices(6, 2) = 1.0;
    vertices(7, 0) = 0.0; vertices(7, 1) = 1.0; vertices(7, 2) = 1.0;
}

template <>
Element<CEL>::Element(TPZFMatrix<REAL> &v)
	: vertices(8, 3),
	coordinates(1, 3),
	idxVec({ 0,1,2,3,4,5 }),
	perm(idxVec)
{
    if (v.Rows() < 8 || v.Cols() < 3)
    {
        throw std::invalid_argument("CubeEl vertices constructor: not enough size in argument (at least 8x3)");
    }
    else
    {
        vertices = v;
    }
}

template <>
Element<CEL>::Element(const Element<CEL> &cp)
	: vertices(cp.vertices),
	coordinates(1, 3),
	idxVec(cp.idxVec),
	perm(cp.perm)
{}

template <>
uint64_t Element<CEL>::position(const TPZVec<uint64_t> &point)
{
    uint64_t n = perm.getPOrder();
    if (point.size() >= 3)
        return perm.getPermutationVector()[ point[0] * (n + 1) * (n + 1) + point[1] * (n + 1) + point[2] ];
    else
        throw new std::logic_error("Cube Element 'position' method called with too few vector elements\n\t3 required");
}

// implemented using nodal shape function
template <>
TPZFMatrix<REAL> &Element<CEL>::mapToElement(TPZFMatrix<REAL> &xi, TPZFMatrix<REAL> &jacobian)
{
    double xi_ = xi(0);
    double eta_ = xi(1);
    double z_ = xi(2);
    {
        double N0 = (1.0 - xi_) * (1.0 - eta_) * (1.0 - z_);
        double N1 = xi_ * (1.0 - eta_) * (1.0 - z_);
        double N2 = xi_ * eta_ * (1.0 - z_);
        double N3 = (1.0 - xi_) * eta_ * (1.0 - z_);
        double N4 = (1.0 - xi_) * (1.0 - eta_) * z_;
        double N5 = xi_ * (1.0 - eta_) * z_;
        double N6 = xi_ * eta_ * z_;
        double N7 = (1.0 - xi_) * eta_ * z_;


        coordinates(0) = N0 * vertices(0, 0) + N1 * vertices(1, 0) + N2 * vertices(2, 0) + N3 * vertices(3, 0) +
            N4 * vertices(4, 0) + N5 * vertices(5, 0) + N6 * vertices(6, 0) + N7 * vertices(7, 0);
        coordinates(1) = N0 * vertices(0, 1) + N1 * vertices(1, 1) + N2 * vertices(2, 1) + N3 * vertices(3, 1) +
            N4 * vertices(4, 1) + N5 * vertices(5, 1) + N6 * vertices(6, 1) + N7 * vertices(7, 1);
        coordinates(2) = N0 * vertices(0, 2) + N1 * vertices(1, 2) + N2 * vertices(2, 2) + N3 * vertices(3, 2) +
            N4 * vertices(4, 2) + N5 * vertices(5, 2) + N6 * vertices(6, 2) + N7 * vertices(7, 2);
    }

    // TODO: implement jacobian
    jacobian(0, 0) = (1.0 - eta_) * (1.0 - z_) * (vertices(1, 0) - vertices(0, 0));
    jacobian(0, 0) += eta_ * (1.0 - z_) * (vertices(2, 0) - vertices(3, 0));
    jacobian(0, 0) += (1.0 - eta_) * z_ * (vertices(5, 0) - vertices(4, 0));
    jacobian(0, 0) += eta_ * z_ * (vertices(6, 0) - vertices(7, 0));

    jacobian(0, 1) = (1.0 - xi_) * (1.0 - z_) * (vertices(3, 0) - vertices(0, 0));
    jacobian(0, 1) += xi_ * (1.0 - z_) * (vertices(2, 0) - vertices(1, 0));
    jacobian(0, 1) += (1.0 - xi_) * z_ * (vertices(7, 0) - vertices(4, 0));
    jacobian(0, 1) += xi_ * z_ * (vertices(6, 0) - vertices(5, 0));

    jacobian(0, 2) = (1.0 - xi_) * (1.0 - eta_) * (vertices(4, 0) - vertices(0, 0));
    jacobian(0, 2) += xi_ * (1.0 - eta_) * (vertices(5, 0) - vertices(1, 0));
    jacobian(0, 2) += (1.0 - xi_) * eta_ * (vertices(6, 0) - vertices(2, 0));
    jacobian(0, 2) += xi_ * eta_ * (vertices(7, 0) - vertices(3, 0));

    jacobian(1, 0) = (1.0 - eta_) * (1.0 - z_) * (vertices(1, 1) - vertices(0, 1));
    jacobian(1, 0) += eta_ * (1.0 - z_) * (vertices(2, 1) - vertices(3, 1));
    jacobian(1, 0) += (1.0 - eta_) * z_ * (vertices(5, 1) - vertices(4, 1));
    jacobian(1, 0) += eta_ * z_ * (vertices(6, 1) - vertices(7, 1));

    jacobian(1, 1) = (1.0 - xi_) * (1.0 - z_) * (vertices(3, 1) - vertices(0, 1));
    jacobian(1, 1) += xi_ * (1.0 - z_) * (vertices(2, 1) - vertices(1, 1));
    jacobian(1, 1) += (1.0 - xi_) * z_ * (vertices(7, 1) - vertices(4, 1));
    jacobian(1, 1) += xi_ * z_ * (vertices(6, 1) - vertices(5, 1));

    jacobian(1, 2) = (1.0 - xi_) * (1.0 - eta_) * (vertices(4, 1) - vertices(0, 1));
    jacobian(1, 2) += xi_ * (1.0 - eta_) * (vertices(5, 1) - vertices(1, 1));
    jacobian(1, 2) += (1.0 - xi_) * eta_ * (vertices(6, 1) - vertices(2, 1));
    jacobian(1, 2) += xi_ * eta_ * (vertices(7, 1) - vertices(3, 1));

    jacobian(2, 0) = (1.0 - eta_) * (1.0 - z_) * (vertices(1, 2) - vertices(0, 2));
    jacobian(2, 0) += eta_ * (1.0 - z_) * (vertices(2, 2) - vertices(3, 2));
    jacobian(2, 0) += (1.0 - eta_) * z_ * (vertices(5, 2) - vertices(4, 2));
    jacobian(2, 0) += eta_ * z_ * (vertices(6, 2) - vertices(7, 2));

    jacobian(2, 1) = (1.0 - xi_) * (1.0 - z_) * (vertices(3, 2) - vertices(0, 2));
    jacobian(2, 1) += xi_ * (1.0 - z_) * (vertices(2, 2) - vertices(1, 2));
    jacobian(2, 1) += (1.0 - xi_) * z_ * (vertices(7, 2) - vertices(4, 2));
    jacobian(2, 1) += xi_ * z_ * (vertices(6, 2) - vertices(5, 2));

    jacobian(2, 2) = (1.0 - xi_) * (1.0 - eta_) * (vertices(4, 2) - vertices(0, 2));
    jacobian(2, 2) += xi_ * (1.0 - eta_) * (vertices(5, 2) - vertices(1, 2));
    jacobian(2, 2) += (1.0 - xi_) * eta_ * (vertices(6, 2) - vertices(2, 2));
    jacobian(2, 2) += xi_ * eta_ * (vertices(7, 2) - vertices(3, 2));

    return coordinates;
}

template <>
TPZFMatrix<REAL>& Element<CEL>::mapToElement(const std::initializer_list<REAL>& xi, TPZFMatrix<REAL>& jacobian)
{
	auto it = xi.begin();

	double xi_ = *it;
	double eta_ = *(++it);
	double z_ = *(++it);
	{
		double N0 = (1.0 - xi_) * (1.0 - eta_) * (1.0 - z_);
		double N1 = xi_ * (1.0 - eta_) * (1.0 - z_);
		double N2 = xi_ * eta_ * (1.0 - z_);
		double N3 = (1.0 - xi_) * eta_ * (1.0 - z_);
		double N4 = (1.0 - xi_) * (1.0 - eta_) * z_;
		double N5 = xi_ * (1.0 - eta_) * z_;
		double N6 = xi_ * eta_ * z_;
		double N7 = (1.0 - xi_) * eta_ * z_;


		coordinates(0) = N0 * vertices(0, 0) + N1 * vertices(1, 0) + N2 * vertices(2, 0) + N3 * vertices(3, 0) +
			N4 * vertices(4, 0) + N5 * vertices(5, 0) + N6 * vertices(6, 0) + N7 * vertices(7, 0);
		coordinates(1) = N0 * vertices(0, 1) + N1 * vertices(1, 1) + N2 * vertices(2, 1) + N3 * vertices(3, 1) +
			N4 * vertices(4, 1) + N5 * vertices(5, 1) + N6 * vertices(6, 1) + N7 * vertices(7, 1);
		coordinates(2) = N0 * vertices(0, 2) + N1 * vertices(1, 2) + N2 * vertices(2, 2) + N3 * vertices(3, 2) +
			N4 * vertices(4, 2) + N5 * vertices(5, 2) + N6 * vertices(6, 2) + N7 * vertices(7, 2);
	}

	// TODO: implement jacobian
	jacobian(0, 0) = (1.0 - eta_) * (1.0 - z_) * (vertices(1, 0) - vertices(0, 0));
	jacobian(0, 0) += eta_ * (1.0 - z_) * (vertices(2, 0) - vertices(3, 0));
	jacobian(0, 0) += (1.0 - eta_) * z_ * (vertices(5, 0) - vertices(4, 0));
	jacobian(0, 0) += eta_ * z_ * (vertices(6, 0) - vertices(7, 0));

	jacobian(0, 1) = (1.0 - xi_) * (1.0 - z_) * (vertices(3, 0) - vertices(0, 0));
	jacobian(0, 1) += xi_ * (1.0 - z_) * (vertices(2, 0) - vertices(1, 0));
	jacobian(0, 1) += (1.0 - xi_) * z_ * (vertices(7, 0) - vertices(4, 0));
	jacobian(0, 1) += xi_ * z_ * (vertices(6, 0) - vertices(5, 0));

	jacobian(0, 2) = (1.0 - xi_) * (1.0 - eta_) * (vertices(4, 0) - vertices(0, 0));
	jacobian(0, 2) += xi_ * (1.0 - eta_) * (vertices(5, 0) - vertices(1, 0));
	jacobian(0, 2) += (1.0 - xi_) * eta_ * (vertices(6, 0) - vertices(2, 0));
	jacobian(0, 2) += xi_ * eta_ * (vertices(7, 0) - vertices(3, 0));

	jacobian(1, 0) = (1.0 - eta_) * (1.0 - z_) * (vertices(1, 1) - vertices(0, 1));
	jacobian(1, 0) += eta_ * (1.0 - z_) * (vertices(2, 1) - vertices(3, 1));
	jacobian(1, 0) += (1.0 - eta_) * z_ * (vertices(5, 1) - vertices(4, 1));
	jacobian(1, 0) += eta_ * z_ * (vertices(6, 1) - vertices(7, 1));

	jacobian(1, 1) = (1.0 - xi_) * (1.0 - z_) * (vertices(3, 1) - vertices(0, 1));
	jacobian(1, 1) += xi_ * (1.0 - z_) * (vertices(2, 1) - vertices(1, 1));
	jacobian(1, 1) += (1.0 - xi_) * z_ * (vertices(7, 1) - vertices(4, 1));
	jacobian(1, 1) += xi_ * z_ * (vertices(6, 1) - vertices(5, 1));

	jacobian(1, 2) = (1.0 - xi_) * (1.0 - eta_) * (vertices(4, 1) - vertices(0, 1));
	jacobian(1, 2) += xi_ * (1.0 - eta_) * (vertices(5, 1) - vertices(1, 1));
	jacobian(1, 2) += (1.0 - xi_) * eta_ * (vertices(6, 1) - vertices(2, 1));
	jacobian(1, 2) += xi_ * eta_ * (vertices(7, 1) - vertices(3, 1));

	jacobian(2, 0) = (1.0 - eta_) * (1.0 - z_) * (vertices(1, 2) - vertices(0, 2));
	jacobian(2, 0) += eta_ * (1.0 - z_) * (vertices(2, 2) - vertices(3, 2));
	jacobian(2, 0) += (1.0 - eta_) * z_ * (vertices(5, 2) - vertices(4, 2));
	jacobian(2, 0) += eta_ * z_ * (vertices(6, 2) - vertices(7, 2));

	jacobian(2, 1) = (1.0 - xi_) * (1.0 - z_) * (vertices(3, 2) - vertices(0, 2));
	jacobian(2, 1) += xi_ * (1.0 - z_) * (vertices(2, 2) - vertices(1, 2));
	jacobian(2, 1) += (1.0 - xi_) * z_ * (vertices(7, 2) - vertices(4, 2));
	jacobian(2, 1) += xi_ * z_ * (vertices(6, 2) - vertices(5, 2));

	jacobian(2, 2) = (1.0 - xi_) * (1.0 - eta_) * (vertices(4, 2) - vertices(0, 2));
	jacobian(2, 2) += xi_ * (1.0 - eta_) * (vertices(5, 2) - vertices(1, 2));
	jacobian(2, 2) += (1.0 - xi_) * eta_ * (vertices(6, 2) - vertices(2, 2));
	jacobian(2, 2) += xi_ * eta_ * (vertices(7, 2) - vertices(3, 2));

	return coordinates;
}