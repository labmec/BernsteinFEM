#include <iostream>
#include <cmath>
#include <cstdlib>
using std::cin;
using std::cout;
using std::endl;

#include "Moments.h"
#include "JacobiGaussNodes.h"

BMoment2DQuad::BMoment2DQuad()
    : Bmoment(), Cval(), quadraWN(), vertices(4, 2, arma::fill::none)
{
    cout << "Enter a value for the polynomial order n:";
    cin >> n;
    cout << endl
         << "Enter a value for the quadrature order q:";
    cin >> q;
    cout << endl;
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }

    quadraWN.set_size(q, 2);
    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);

    Bmoment.zeros(lenMoments, 0);
    Cval.set_size(q * q, 1);
}

BMoment2DQuad::BMoment2DQuad(int q, int n)
    : Bmoment((MAX(n + 1, q) * MAX(n + 1, q)), 1, arma::fill::zeros),
      Cval(q * q, 1, arma::fill::none),
      quadraWN(q, 2, arma::fill::none),
      vertices(4, 2, arma::fill::none)
{
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;
    this->m = n;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);
}

BMoment2DQuad::BMoment2DQuad(int q, int n, int m, int nb_Array)
    : Bmoment(MAX(n + 1, m + 1) * MAX(n + 1, m + 1), nb_Array, arma::fill::zeros),
      Cval(q * q, nb_Array, arma::fill::none),
      quadraWN(q, 2, arma::fill::none),
      vertices(4, 2, arma::fill::none)
{
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;
    this->m = m;

    assignQuadra();

    int mx = MAX(n + 1, m + 1);
    lenMoments = mx * mx;
}

BMoment2DQuad::~BMoment2DQuad()
{
}

void BMoment2DQuad::assignQuadra()
{
    double *x = quadraWN.colptr(1);
    double *w = quadraWN.colptr(0);

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

void BMoment2DQuad::nodalShape(double X[], double &dX, double xi, double eta)
{
    arma::vec x(2, arma::fill::none);
    nodalShape(x, dX, xi, eta);
}

void BMoment2DQuad::nodalShape(arma::vec &X, double &dX, double xi, double eta)
{
    arma::mat a(2, 2, arma::fill::none);
    nodalShape(X, a, dX, xi, eta);
}

void BMoment2DQuad::nodalShape(arma::vec &X, arma::mat &jac, double &dX, double xi, double eta)
{
    double N[4];
    double x_xi, x_eta, y_xi, y_eta;

    //computes the nodal-shape values
    N[0] = (1.0 - xi) * (1.0 - eta);
    N[1] = xi * (1.0 - eta);
    N[2] = xi * eta;
    N[3] = (1.0 - xi) * eta;

    //computes the mapped point
    X(0) = N[0] * vertices(0, 0) + N[1] * vertices(1, 0) + N[2] * vertices(2, 0) + N[3] * vertices(3, 0);
    X(1) = N[0] * vertices(0, 1) + N[1] * vertices(1, 1) + N[2] * vertices(2, 1) + N[3] * vertices(3, 1);

    //computes the derivatives
    x_xi = (1.0 - eta) * (vertices(1, 0) - vertices(0, 0)) + eta * (vertices(2, 0) - vertices(3, 0));
    x_eta = (1.0 - xi) * (vertices(3, 0) - vertices(0, 0)) + xi * (vertices(2, 0) - vertices(1, 0));
    y_xi = (1.0 - eta) * (vertices(1, 1) - vertices(0, 1)) + eta * (vertices(2, 1) - vertices(3, 1));
    y_eta = (1.0 - xi) * (vertices(3, 1) - vertices(0, 1)) + xi * (vertices(2, 1) - vertices(1, 1));

    // stores the Jacobian matrix
    jac(0, 0) = x_xi;
    jac(0, 1) = x_eta;
    jac(1, 0) = y_xi;
    jac(1, 1) = y_eta;

    //computes the Jacobian det
    dX = x_xi * y_eta - x_eta * y_xi;
}

void BMoment2DQuad::computeFunctionDef()
{
    int i, j;
    double X[2], dX;
    int index_ij;

    for (i = 0; i < q; i++)
    {
        for (j = 0; j < q; j++)
        {
            nodalShape(X, dX, quadraWN(i, 1), quadraWN(j, 1));
            index_ij = position(i, j, q);
            Cval(index_ij, 0) = (*f)(X[0], X[1]) * dX;
        }
    }

    fValSet = true;
}

// it is expected that Fval has the function values mapped to the
// master element already multiplied by the Jacobian
// and has q * nb_Array elements
void BMoment2DQuad::setFunction(arma::vec Fval)
{
    for (int i = 0; i < q * q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval(i, el) = Fval(i + el * q);

    fValSet = true;
}
// Fval must have at least q X nb_Array elements
void BMoment2DQuad::setFunction(arma::mat Fval)
{
    for (int i = 0; i < q * q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval(i, el) = Fval(i, el);

    fValSet = true;
}

void BMoment2DQuad::setFunction(double (*function)(double, double))
{
    f = function;
    fDefSet = true;
}

void BMoment2DQuad::setQuadrilateral(double v1[2], double v2[2], double v3[2], double v4[2])
{
    vertices(0, 0) = v1[0];
    vertices(0, 1) = v1[1];
    vertices(1, 0) = v2[0];
    vertices(1, 1) = v2[1];
    vertices(2, 0) = v3[0];
    vertices(2, 1) = v3[1];
    vertices(3, 0) = v4[0];
    vertices(3, 1) = v4[1];
}

void BMoment2DQuad::setQuadrilateral(arma::vec v1, arma::vec v2, arma::vec v3, arma::vec v4)
{
    vertices(0, 0) = v1(0);
    vertices(0, 1) = v1(1);
    vertices(1, 0) = v2(0);
    vertices(1, 1) = v2(1);
    vertices(2, 0) = v3(0);
    vertices(2, 1) = v3(1);
    vertices(3, 0) = v4(0);
    vertices(3, 1) = v4(1);
}

void BMoment2DQuad::setQuadrilateral(arma::mat vertices)
{
    this->vertices = vertices;
}

arma::mat BMoment2DQuad::getIntegrationPoints()
{
    arma::mat points(q * q, 2);
    arma::vec X(2, arma::fill::none);
    double dX;

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            nodalShape(X, dX, quadraWN(i, 1), quadraWN(j, 1));
            points(i * q + j, 0) = X[0];
            points(i * q + j, 1) = X[1];
        }
    }

    return points;
}

void BMoment2DQuad::compute_moments()
{
    if (functVal == 1 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 0 && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        if (functVal == 0)
            computeFunctionDef();

        double xi1, xi2, omega1, omega2, s1, s2, r1, r2, w1, w2;
        int i, j, a1, a2, index_ij, index_a1a2;

        for (i = 0; i < q; i++)
        {
            xi1 = quadraWN(i, 1);    //quadrature abcissa i
            omega1 = quadraWN(i, 0); //quadrature weight i
            s1 = 1.0 - xi1;
            r1 = xi1 / s1; //recurrence relation 1st coefficient
            for (a1 = 0; a1 <= n; a1++)
            {
                w1 = omega1 * pow(s1, n);
                for (j = 0; j < q; j++)
                {
                    xi2 = quadraWN(j, 1);    //quadrature abcissa j
                    omega2 = quadraWN(j, 0); //quadrature weight j
                    s2 = 1.0 - xi2;
                    r2 = xi2 / s2; //recurrence relation 2nd coefficient

                    index_ij = position(i, j, q - 1);

                    w2 = omega2 * pow(s2, n);
                    for (a2 = 0; a2 <= m; a2++)
                    {
                        // here w2 equals to the weight j multiplied with the a2-th 1d bernstein polynomial evaluated at xi2
                        index_a1a2 = position(a1, a2, m);
                        for (int el = 0; el < nb_Array; el++)
                            Bmoment(index_a1a2, el) += w1 * w2 * Cval(index_ij, el);
                        w2 *= r2 * ((m - a2) / (1.0 + a2));
                    }
                }
                w1 *= r1 * ((n - a1) / (1.0 + a1));
            }
        }
    }
}

void BMoment2DQuad::compute_moments(double (*f)(double, double))
{
    setFunction(f);
    compute_moments();
}

void BMoment2DQuad::compute_moments(arma::vec Fval)
{
    setFunction(Fval);
    compute_moments();
}