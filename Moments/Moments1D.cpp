#include <iostream>
#include <cmath>
#include <cstdlib>
using std::cin;
using std::cout;
using std::endl;

#include "Moments.h"
#include "JacobiGaussNodes.h"

// default constructor
BMoment1D::BMoment1D()
    : Bmoment(), Cval(), quadraWN()
{
    cout << "Enter a value for the polynomial order n:";
    cin >> n;
    cout << endl
         << "Enter a value for the quadrature order q:";
    cin >> q;
    cout << endl;
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
    }

    int m = MAX(n, q - 1);
    lenMoments = (m + 1);

    Bmoment.zeros(lenMoments, 1);
    Cval.zeros(q, 1);
}

// quadrature and polynomial order constructor
BMoment1D::BMoment1D(int q, int n)
    : Bmoment(MAX(n + 1, q), 1, arma::fill::zeros), Cval(q, 1, arma::fill::zeros),
      quadraWN(q + 1, 2, arma::fill::zeros)
{
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The number of quadrature points is too large.\n";
        throw std::bad_alloc();
    }
    this->q = q;
    this->n = n;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1);
}

BMoment1D::BMoment1D(int q, int n, int nb_Array)
    : Bmoment(MAX(n + 1, q), nb_Array, arma::fill::zeros),
      Cval(q, nb_Array, arma::fill::zeros),
      quadraWN(q + 1, 2, arma::fill::zeros)
{
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The number of quadrature points is too large.\n";
        throw std::bad_alloc();
    }
    this->q = q;
    this->n = n;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1);
}

BMoment1D::~BMoment1D()
{
}

void BMoment1D::assignQuadra()
{
    double *x = quadraWN.colptr(1);
    double *w = quadraWN.colptr(0);

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

void BMoment1D::loadFunctionDef()
{
    for (int i = 0; i < q; i++)
    {
        for (int el = 0; el < nb_Array; el++)
            Cval.at(i, el) = f(quadraWN(i, 0));
    }

    fValSet = true;
}

// set the function values for computation (Fval must have at least q * nb_Array elements)
void BMoment1D::setFunction(const arma::vec &Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval.at(i, el) = Fval.at(i + el * q);

    fValSet = true;
}

// Fval must have at least q X nb_Array elements
void BMoment1D::setFunction(const arma::mat &Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval.at(i, el) = Fval.at(i, el);

    fValSet = true;
}

// set the function definition for computation
void BMoment1D::setFunction(std::function<double (double)> function)
{
    f = function;
    fDefSet = true;
}

void BMoment1D::setInterval(double a, double b)
{
    this->a = a;
    this->b = b;
}

// returns the vector with the integration points over the object's element, following the moments organization
arma::vec BMoment1D::getIntegrationPoints()
{
    arma::vec points(this->q); // vector with q elements

    for(int i = 0; i < q; i++)
    {
        points.at(i) = (legendre_xi(q, i) + (a + 1)) * ((b - a) / 2.);
    }

    return points;
}

// compute the B-moments using the values already assigned in the object
void BMoment1D::compute_moments()
{
    if (functVal == 0 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 1 && !fValSet)
    {
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
        loadFunctionDef();
    }
    else
    {
        Bmoment.zeros();
        for (int i = 0; i < q; i++)
        {
            double xi = quadraWN.at(i, 1);
            double omega = quadraWN.at(i, 0);
            double s = 1 - xi;
            double r = xi / s;
            double w = omega * pow(s, n);
            for (int alpha = 0; alpha <= n; alpha++)
            {
                // here w equals the Bernstein polynomial of order n,
                // with index alpha, evaluated at the i-th integration node
                // times the i-th integration weight.
                int p = position(alpha, n);
                for (int el = 0; el < nb_Array; el++)
                    Bmoment.at(p, el) += w * Cval.at(i, el);
                w *= r * ((n - alpha) / (1. + alpha)); // treats the recurrence relation
            }
        }
    }
}

// compute the b-moments for the specified f function
void BMoment1D::compute_moments(std::function<double (double)> f)
{
    setFunction(f);
    compute_moments();
}

// compute the b-moments for the Fval function values
void BMoment1D::compute_moments(const arma::vec &Fval)
{
    setFunction(Fval);
    compute_moments();
}
