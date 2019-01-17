#include "Elements.h"

#define TEL Element_t::TetrahedronEl

template <>
arma::mat Element<TEL>::jac(3, 3, arma::fill::zeros);

template <>
Element<TEL>::Element()
    : vertices(4, 3, arma::fill::none),
      coordinates(1, 3, arma::fill::none)
{
    vertices.at(0, 0) = 0.0; vertices.at(0, 1) = 0.0; vertices.at(0, 2) = 0.0; // v1 = (0,0,0)
    vertices.at(1, 0) = 1.0; vertices.at(1, 1) = 0.0; vertices.at(1, 2) = 0.0; // v2 = (1,0,0)
    vertices.at(2, 0) = 0.0; vertices.at(2, 1) = 1.0; vertices.at(2, 2) = 0.0; // v3 = (0,1,0)
    vertices.at(3, 0) = 0.0; vertices.at(3, 1) = 0.0; vertices.at(3, 2) = 1.0; // v4 = (0,0,1)
}

template <>
Element<TEL>::Element(const arma::mat &v)
    : vertices(3, 2, arma::fill::none),
      coordinates(1, 2, arma::fill::none)
{
    if (v.n_rows < 3 || v.n_cols < 2)
    {
        throw std::invalid_argument("TetrahedronEl vertices constructor: not enough size in matrix (at least 3x2)");
    }
    else
    {
        vertices.at(0, 0) = v(0, 0); vertices.at(0, 1) = v(0, 1); vertices.at(0, 2) = v(0, 2);
        vertices.at(1, 0) = v(1, 0); vertices.at(1, 1) = v(1, 1); vertices.at(1, 2) = v(1, 2);
        vertices.at(2, 0) = v(2, 0); vertices.at(2, 1) = v(2, 1); vertices.at(2, 2) = v(2, 2);
        vertices.at(3, 0) = v(3, 0); vertices.at(3, 1) = v(3, 1); vertices.at(3, 2) = v(3, 2);
    }
}

template <>
Element<TEL>::Element(const Element<TEL> &cp)
    : vertices(cp.vertices),
      coordinates(1, 2, arma::fill::none) {}

// Duffy Transform
template <>
const arma::mat &Element<TEL>::mapToElement(const arma::mat &xi, arma::mat &jacobian)
{
    try 
    {
        // treat these as barycentric coordinates
        double b1 = xi(0);
        double b2 = xi(1) * (1 - b1);
        double b3 = xi(2) * (1 - b2) * (1 - b1);
        double b4 = 1 - b3 - b2 - b1;
        
        // convert from barycentric coordinates to cartesian
        coordinates(0) = b1 * vertices.at(0, 0) + b2 * vertices.at(1, 0) +
            b3 * vertices.at(2, 0) + b4 * vertices.at(3, 0);
        coordinates(1) = b1 * vertices.at(0, 1) + b2 * vertices.at(1, 1) +
            b3 * vertices.at(2, 1) + b4 * vertices.at(3, 1);
        coordinates(2) = b1 * vertices.at(0, 2) + b2 * vertices.at(1, 2) +
            b3 * vertices.at(2, 2) + b4 * vertices.at(3, 2);
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