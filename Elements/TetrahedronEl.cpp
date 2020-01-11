#include "Elements.h"

#define TEL Element_t::TetrahedronEl

template <>
TPZFMatrix<REAL> Element<TEL>::jac(3, 3, 0);

template <>
Element<TEL>::Element()
	: vertices(4, 3),
	coordinates(1, 3),
	idxVec({ 0,1,2,3 }),
	perm(idxVec)
{
    vertices(0, 0) = 0.0; vertices(0, 1) = 0.0; vertices(0, 2) = 0.0; // v1 = (0,0,0)
    vertices(1, 0) = 1.0; vertices(1, 1) = 0.0; vertices(1, 2) = 0.0; // v2 = (1,0,0)
    vertices(2, 0) = 0.0; vertices(2, 1) = 1.0; vertices(2, 2) = 0.0; // v3 = (0,1,0)
    vertices(3, 0) = 0.0; vertices(3, 1) = 0.0; vertices(3, 2) = 1.0; // v4 = (0,0,1)
}

template <>
Element<TEL>::Element(TPZFMatrix<REAL> &v)
	: vertices(3, 2),
	coordinates(1, 2),
	idxVec({ 0,1,2,3 }),
	perm(idxVec)
{
    if (v.Rows() < 3 || v.Cols() < 2)
    {
        throw std::invalid_argument("TetrahedronEl vertices constructor: not enough size in matrix (at least 3x2)");
    }
    else
    {
        vertices(0, 0) = v(0, 0); vertices(0, 1) = v(0, 1); vertices(0, 2) = v(0, 2);
        vertices(1, 0) = v(1, 0); vertices(1, 1) = v(1, 1); vertices(1, 2) = v(1, 2);
        vertices(2, 0) = v(2, 0); vertices(2, 1) = v(2, 1); vertices(2, 2) = v(2, 2);
        vertices(3, 0) = v(3, 0); vertices(3, 1) = v(3, 1); vertices(3, 2) = v(3, 2);
    }
}

template <>
uint64_t Element<TEL>::position(const TPZVec<uint64_t> &point)
{
    const uint64_t n = perm.getPOrder();
    if (point.size() >= 3)
        return perm.getPermutationVector()[ point[0] * (n + 1) * (n + 1) + point[1] * (n + 1) + point[2] ];
    else
        throw std::logic_error("Tetrahedron Element 'position' method called with too few vector elements\n\t3 required");
    ;
}

template <>
uint64_t Element<TEL>::position(const std::initializer_list<uint64_t> &point)
{
    const uint64_t n = perm.getPOrder();
    if (point.size() >= 3) {
	    auto it = point.begin();
    	const uint64_t a1 = *it, a2 = *(++it), a3 = *(++it);
        return perm.getPermutationVector()[ a1 * (n + 1) * (n + 1) + a2 * (n + 1) + a3 ];
    } else {
        throw std::logic_error("Tetrahedron Element 'position' method called with too few vector elements\n\t3 required");
    }
}

template <>
Element<TEL>::Element(const Element<TEL> &cp)
	: vertices(cp.vertices),
	coordinates(1, 2),
	idxVec(cp.idxVec),
	perm(cp.perm)
{}

// Duffy Transform
template <>
TPZFMatrix<REAL> &Element<TEL>::mapToElement(TPZFMatrix<REAL> &xi, TPZFMatrix<REAL> &jacobian)
{
    try 
    {
        // treat these as barycentric coordinates
        double b1 = xi(0);
        double b2 = xi(1) * (1 - b1);
        double b3 = xi(2) * (1 - b2) * (1 - b1);
        double b4 = 1 - b3 - b2 - b1;
        
        // convert from barycentric coordinates to cartesian
        coordinates(0) = b1 * vertices(0, 0) + b2 * vertices(1, 0) +
            b3 * vertices(2, 0) + b4 * vertices(3, 0);
        coordinates(1) = b1 * vertices(0, 1) + b2 * vertices(1, 1) +
            b3 * vertices(2, 1) + b4 * vertices(3, 1);
        coordinates(2) = b1 * vertices(0, 2) + b2 * vertices(1, 2) +
            b3 * vertices(2, 2) + b4 * vertices(3, 2);
    }
    catch(std::logic_error &e)
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