#include "Elements.h"

#define TEL Element_t::TriangularEl

template <>
TPZFMatrix<REAL> Element<TEL>::jac(2, 2, 0);

template <>
Element<TEL>::Element()
	: vertices(3, 2),
	coordinates(1, 2),
	idxVec({ 0,1,2 }),
	perm(idxVec)
{
    vertices(0, 0) = 0.0; vertices(0, 1) = 0.0; // v1 = (0,0)
    vertices(1, 0) = 1.0; vertices(1, 1) = 0.0; // v2 = (1,0)
    vertices(2, 0) = 0.0; vertices(2, 1) = 1.0; // v3 = (0,1)
}

template <>
Element<TEL>::Element(TPZFMatrix<REAL> &v)
	: vertices(3, 2),
	coordinates(1, 2),
	idxVec({ 0,1,2 }),
	perm(idxVec)
{
    if (v.Rows() < 3 || v.Cols() < 2)
    {
        throw std::invalid_argument("TriangularEl vertices constructor: not enough size in matrix (at least 3x2)");
    }
    else
    {
        vertices(0, 0) = v(0, 0); vertices(0, 1) = v(0, 1);
        vertices(1, 0) = v(1, 0); vertices(1, 1) = v(1, 1);
        vertices(2, 0) = v(2, 0); vertices(2, 1) = v(2, 1);
    }
}

template <>
Element<TEL>::Element(const Element<TEL> &cp)
	: vertices(cp.vertices),
	coordinates(1, 2),
	idxVec(cp.idxVec),
	perm(cp.perm) 
{}

template <>
uint64_t Element<TEL>::position(const TPZVec<uint64_t> &point)
{
    uint64_t n = perm.getPOrder();
    if (point.size() >= 2)
        return perm.getPermutationVector()[ point[0] * (n + 1) + point[1] ];
    else
        throw std::logic_error("Triangular Element 'position' method called with too few vector elements\n\t2 required");
}

template <>
uint64_t Element<TEL>::position(const std::initializer_list<uint64_t>& point)
{
	uint64_t n = perm.getPOrder();
	if (point.size() >= 2)
	{
		auto it = point.begin();
		uint64_t a1 = *it, a2 = *(++it);
		return perm.getPermutationVector()[a1 * (n + 1) + a2];
	} else {
		throw std::logic_error("Triangular Element 'position' method called with too few vector elements\n\t2 required");
	}
}

// Duffy Transform
template <>
TPZFMatrix<REAL> &Element<TEL>::mapToElement(TPZFMatrix<REAL> &xi, TPZFMatrix<REAL> &jacobian)
{
    try 
    {
        // treat these as barycentric coordinates
        double b1 = xi(0);
        double b2 = xi(1) * (1 - b1);
        double b3 = 1 - b2 - b1;
        
        // convert from barycentric coordinates to cartesian
        coordinates(0) = b1 * vertices(0, 0) + b2 * vertices(1, 0) + b3 * vertices(2, 0);
        coordinates(1) = b1 * vertices(0, 1) + b2 * vertices(1, 1) + b3 * vertices(2, 1);
    }
    catch(std::logic_error)
    {
        throw std::invalid_argument("TriangularEl mapToElement first argument: not enough size (at least 2)");
    }
    catch (std::exception &e)
    {
        std::cerr << e.what() << ", in method Element<TriangularEl>::mapToElement" << std::endl;
        std::terminate();
    }

    // TODO: implement jacobian

    return coordinates;
}

template <>
TPZFMatrix<REAL>& Element<TEL>::mapToElement(const std::initializer_list<REAL> &xi, TPZFMatrix<REAL> &jacobian)
{
	try
	{
		auto it = xi.begin();
		// treat these as barycentric coordinates
		double b1 = *it;
		double b2 = *(++it) * (1 - b1);
		double b3 = 1 - b2 - b1;

		// convert from barycentric coordinates to cartesian
		coordinates(0) = b1 * vertices(0, 0) + b2 * vertices(1, 0) + b3 * vertices(2, 0);
		coordinates(1) = b1 * vertices(0, 1) + b2 * vertices(1, 1) + b3 * vertices(2, 1);
	}
	catch (std::logic_error)
	{
		throw std::invalid_argument("TriangularEl mapToElement first argument: not enough size (at least 2)");
	}
	catch (std::exception& e)
	{
		std::cerr << e.what() << ", in method Element<TriangularEl>::mapToElement" << std::endl;
		std::terminate();
	}

	return coordinates;
}