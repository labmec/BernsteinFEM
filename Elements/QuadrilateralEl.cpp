#include "Elements.h"

#define QEL Element_t::QuadrilateralEl

template <>
arma::mat Element<QEL>::jac(2, 2, arma::fill::zeros);

template <>
Element<QEL>::Element()
    : vertices(4, 2, arma::fill::none),
      coordinates(1, 2, arma::fill::none)
{
    vertices(0, 0) = 0.0; vertices(0, 1) = 0.0; // v1 = (0,0)
    vertices(1, 0) = 1.0; vertices(1, 1) = 0.0; // v2 = (1,0)
    vertices(2, 0) = 1.0; vertices(2, 1) = 1.0; // v3 = (1,1)
    vertices(3, 0) = 0.0; vertices(3, 1) = 1.0; // v4 = (0,1)
}

template <>
Element<QEL>::Element(const arma::mat &v)
    : vertices(4, 2, arma::fill::none),
      coordinates(1, 2, arma::fill::none)
{
    if (v.n_rows < 4 || v.n_cols < 2)
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
      coordinates(1, 2, arma::fill::none) {}

template <>
uint Element<QEL>::position(const std::vector<uint> &point, int n)
{
    if (point.size() >= 2)
        return point[0] * (n + 1) + point[1];
    else
        throw new std::logic_error("Quadrilateral Element 'position' method called with too few vector elements\n\t2 required");
}

// implemented using nodal shape function
template <>
const arma::mat &Element<QEL>::mapToElement(const arma::mat &xi, arma::mat &jacobian)
{
    double xi_ = xi(0);
    double eta_ = xi(1);
    {
        double N0 = (1.0 - xi_) * (1.0 - eta_);
        double N1 = xi_ * (1.0 - eta_);
        double N2 = xi_ * eta_;
        double N3 = (1.0 - xi_) * eta_;

        coordinates(0) = N0 * vertices.at(0, 0) + N1 * vertices.at(1, 0) + N2 * vertices.at(2, 0) + N3 * vertices.at(3, 0);
        coordinates(1) = N0 * vertices.at(0, 1) + N1 * vertices.at(1, 1) + N2 * vertices.at(2, 1) + N3 * vertices.at(3, 1);
    }

    jacobian.at(0, 0) = (1.0 - eta_) * (vertices.at(1, 0) - vertices.at(0, 0)) + eta_ * (vertices.at(2, 0) - vertices.at(3, 0));
    jacobian.at(0, 1) = (1.0 - xi_) * (vertices.at(3, 0) - vertices.at(0, 0)) + xi_ * (vertices.at(2, 0) - vertices.at(1, 0));
    jacobian.at(1, 0) = (1.0 - eta_) * (vertices.at(1, 1) - vertices.at(0, 1)) + eta_ * (vertices.at(2, 1) - vertices.at(3, 1));
    jacobian.at(1, 1) = (1.0 - xi_) * (vertices.at(3, 1) - vertices.at(0, 1)) + xi_ * (vertices.at(2, 1) - vertices.at(1, 1));

    return coordinates;
}