#include "Elements.h"

template <>
arma::mat Element<LinearEl>::jac(1, 1, arma::fill::zeros);

template <>
Element<LinearEl>::Element()
    : vertices(2, 1, arma::fill::none),
      coordinates(1, 1, arma::fill::none)
{
    vertices(0) = 0.0;
    vertices(1) = 1.0;
}

template <>
Element<LinearEl>::Element(const arma::mat &v)
    : vertices(2, 1, arma::fill::none),
      coordinates(1, 1, arma::fill::none)
{
    if (vertices.n_rows < 2)
    {
        throw std::invalid_argument("Invalid constructor argument, LinearEl vertices should have at least size 2");
    }
    else
    {
        vertices(0) = v(0);
        vertices(1) = v(1);
    }
}

template <>
Element<LinearEl>::Element(const Element<LinearEl> &cp)
    : vertices(cp.vertices),
      coordinates(1, 1, arma::fill::none) {}

template <>
const arma::mat &Element<LinearEl>::mapToElement(const arma::mat &coordinates, arma::mat &jacobian)
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
}