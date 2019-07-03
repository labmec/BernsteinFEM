#include "Elements.h"

#define LEL Element_t::LinearEl

template <>
TPZFMatrix<REAL> Element<LEL>::jac(1, 1, 0);

template <>
Element<LEL>::Element()
    : vertices(2),
      coordinates(1),
      idxVec({0, 1}),
      perm(idxVec)
{
    vertices(0) = 0.0;
    vertices(1) = 1.0;
}

template <>
Element<LEL>::Element(const TPZFMatrix<REAL> &v)
    : vertices(2),
      coordinates(1),
      idxVec({0, 1}),
      perm(idxVec)
{
    if (vertices.Rows() < 2)
    {
        throw std::invalid_argument("LinearEl vertices constructor: not enough size in argument (at least 2)");
    }
    else
    {
        vertices(0) = v(0);
        vertices(1) = v(1);
    }
}

template <>
Element<LEL>::Element(const Element<LEL> &cp)
    : vertices(cp.vertices),
      coordinates(1),
      idxVec(cp.idxVec),
      perm(cp.perm)
{
}

template <>
uint Element<LEL>::position(const std::vector<uint> &point)
{
    if (point.size() >= 1)
        return perm.getPermutationVector()[point[0]];
    else
        throw new std::logic_error("Linear Element 'position' method called with too few vector elements\n\t1 required");
}

template <>
const TPZFMatrix<REAL> &Element<LEL>::mapToElement(const TPZFMatrix<REAL> &coordinates, TPZFMatrix<REAL> &jacobian)
{
    try
    {
        this->coordinates(0) = (coordinates(0) + vertices(0)) * (vertices(1) - vertices(0));
        jacobian(0) = vertices(1) - vertices(0);
    }
    catch (std::logic_error &e)
    {
        std::cerr << e.what() << std::endl;
        std::cerr << "invalid arguments (probably the size) in 'mapToElement' method call" << std::endl;
    }

    return coordinates;
}