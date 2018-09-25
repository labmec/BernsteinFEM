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
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
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
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
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
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
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
    X(0) = N[0] * vertices.at(0, 0) + N[1] * vertices.at(1, 0) + N[2] * vertices.at(2, 0) + N[3] * vertices.at(3, 0);
    X(1) = N[0] * vertices.at(0, 1) + N[1] * vertices.at(1, 1) + N[2] * vertices.at(2, 1) + N[3] * vertices.at(3, 1);

    //computes the derivatives
    x_xi = (1.0 - eta) * (vertices.at(1, 0) - vertices.at(0, 0)) + eta * (vertices.at(2, 0) - vertices.at(3, 0));
    x_eta = (1.0 - xi) * (vertices.at(3, 0) - vertices.at(0, 0)) + xi * (vertices.at(2, 0) - vertices.at(1, 0));
    y_xi = (1.0 - eta) * (vertices.at(1, 1) - vertices.at(0, 1)) + eta * (vertices.at(2, 1) - vertices.at(3, 1));
    y_eta = (1.0 - xi) * (vertices.at(3, 1) - vertices.at(0, 1)) + xi * (vertices.at(2, 1) - vertices.at(1, 1));

    // stores the Jacobian matrix
    jac.at(0, 0) = x_xi;
    jac.at(0, 1) = x_eta;
    jac.at(1, 0) = y_xi;
    jac.at(1, 1) = y_eta;

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
            nodalShape(X, dX, quadraWN.at(i, 1), quadraWN.at(j, 1));
            index_ij = position(i, j, q);
            Cval.at(index_ij, 0) = f(X[0], X[1]) * dX;
        }
    }

    fValSet = true;
}

// it is expected that Fval has the function values mapped to the
// master element already multiplied by the Jacobian
// and has q * nb_Array elements
void BMoment2DQuad::setFunction(const arma::vec &Fval)
{
    for (int i = 0; i < q * q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval.at(i, el) = Fval.at(i + el * q);

    fValSet = true;
}
// Fval must have at least q X nb_Array elements
void BMoment2DQuad::setFunction(const arma::mat &Fval)
{
    for (int i = 0; i < q * q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval.at(i, el) = Fval.at(i, el);

    fValSet = true;
}

void BMoment2DQuad::setFunction(std::function<double (double, double)> function)
{
    f = function;
    fDefSet = true;
}

void BMoment2DQuad::setQuadrilateral(double v1[2], double v2[2], double v3[2], double v4[2])
{
    vertices.at(0, 0) = v1[0];
    vertices.at(0, 1) = v1[1];
    vertices.at(1, 0) = v2[0];
    vertices.at(1, 1) = v2[1];
    vertices.at(2, 0) = v3[0];
    vertices.at(2, 1) = v3[1];
    vertices.at(3, 0) = v4[0];
    vertices.at(3, 1) = v4[1];
}

void BMoment2DQuad::setQuadrilateral(const arma::vec &v1, const arma::vec &v2, const arma::vec &v3, const arma::vec &v4)
{
    vertices.at(0, 0) = v1.at(0);
    vertices.at(0, 1) = v1.at(1);
    vertices.at(1, 0) = v2.at(0);
    vertices.at(1, 1) = v2.at(1);
    vertices.at(2, 0) = v3.at(0);
    vertices.at(2, 1) = v3.at(1);
    vertices.at(3, 0) = v4.at(0);
    vertices.at(3, 1) = v4.at(1);
}

void BMoment2DQuad::setQuadrilateral(const arma::mat &vertices)
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
            nodalShape(X, dX, quadraWN.at(i, 1), quadraWN.at(j, 1));
            points.at(i * q + j, 0) = X[0];
            points.at(i * q + j, 1) = X[1];
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

        int max_nq = MAX(n, q - 1);
        int max_nm = MAX(n, m);
        arma::mat Bmoment_inter((max_nq + 1) * (max_nq + 1), 1, arma::fill::zeros);
        Bmoment.zeros();

        // convert first index (l=2)
        for (int i = 0; i < q; i++)
        {
            double xi = quadraWN(i, 1);
            double wgt = quadraWN(i, 0);

            double s = 1 - xi;
            double r = xi / s;

            double B = wgt * pow(s, n);
            for (int a1 = 0; a1 <= n; a1++)
            {
                for (int j = 0; j < q; j++)
                {
                    int index_a1j = position(a1, j, max_nq);

                    int index_ij = position(i, j, q - 1);

                    for (int ell = 0; ell < nb_Array; ell++)
                    try {
                        Bmoment_inter(index_a1j, ell) += B * Cval(index_ij, ell);
                    } catch (std::exception e) {
                        std::cout << e.what() << std::endl;
                        std::cout << "index_a1j = " << index_a1j << "MAX: " << Bmoment_inter.n_rows << std::endl;
                        std::cout << "index_ij = " << index_ij << "MAX: " << Cval.n_rows << std::endl;
                        std::terminate();
                    }
                    
                }
                B = B * r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=1)
        for (int i = 0; i < q; i++)
        {
            double xi = quadraWN(i, 1);
            double wgt = quadraWN(i, 0);

            double s = 1 - xi;
            double r = xi / s;

            for (int a1 = 0; a1 <= n; a1++)
            {
                double B = wgt * pow(s, m);
                for (int a2 = 0; a2 <= m; a2++)
                {
                    int index_a1a2 = position(a1, a2, max_nm);
                    int index_a1i = position(a1, i, max_nq);

                    for (int ell = 0; ell < nb_Array; ell++)
                    {
                        try {
                            Bmoment(index_a1a2, 0) += B * Bmoment_inter(index_a1i, 0);
                        } catch (std::exception e ) {
                            std::cout << e.what() << std::endl;
                            std::cout << "index_a1a2 = " << index_a1a2 << "MAX: " << Bmoment.n_rows << std::endl;
                            std::cout << "index_a1i = " << index_a1i << "MAX: " << Bmoment_inter.n_rows << std::endl;
                            std::terminate();
                        }
                    }
                    
                    B = B * r * (m - a2) / (1 + a2);
                }
            }
        }

    }
}

void BMoment2DQuad::compute_moments(std::function<double (double, double)> f)
{
    setFunction(f);
    compute_moments();
}

void BMoment2DQuad::compute_moments(const arma::vec &Fval)
{
    setFunction(Fval);
    compute_moments();
}