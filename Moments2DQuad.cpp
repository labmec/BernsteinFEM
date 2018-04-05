#include <iostream>
#include <cmath>
#include <cstdlib>
using std::cin;
using std::cout;
using std::endl;

#include "Moments.h"
#include "JacobiGaussNodes.h"

BMoment2DQuad::BMoment2DQuad()
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

    quadraWN = create_quadraWN();
    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);

    Bmoment = create_Bmoment();
    Cval = create_Cval();
}

BMoment2DQuad::BMoment2DQuad(int q, int n)
{
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;

    quadraWN = create_quadraWN();
    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);

    Bmoment = create_Bmoment();
    Cval = create_Cval();
}

BMoment2DQuad::~BMoment2DQuad()
{
    delete_Bmoment(Bmoment);
    delete quadraWN[0];
    delete quadraWN;
    delete Cval[0];
    delete Cval;
}

double **BMoment2DQuad::create_Bmoment()
{
    double **Bmoment = new double *[lenMoments];
    double *Bmoment_array = new double[nb_Array * lenMoments];

    double *p1 = Bmoment_array;

    for (int i = 0; i < lenMoments; p1 += nb_Array, i++)
        Bmoment[i] = p1;

    return Bmoment;
}

double **BMoment2DQuad::create_Cval()
{
    int aux = lenMoments;
    lenMoments = q;
    double **p = create_Bmoment();
    lenMoments = aux;
    return p;
}

void BMoment2DQuad::delete_Bmoment(double **Bmoment)
{
    delete Bmoment[0];
    delete Bmoment;
}

double **BMoment2DQuad::create_quadraWN()
{
    double *aux = new double[2 * (q + 1)];
    double **quadraWN = new double *[2];
    quadraWN[0] = aux;
    quadraWN[1] = aux + q + 1;
    return quadraWN;
}

void BMoment2DQuad::assignQuadra()
{
    double *x = quadraWN[1];
    double *w = quadraWN[0];

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

void BMoment2DQuad::nodalShape(double X[], double &dX, double xi, double eta)
{
    double N[4];
    double x_xi, x_eta, y_xi, y_eta;

    //computes the nodal-shape values
    N[0] = (1.0 - xi) * (1.0 - eta);
    N[1] = xi * (1.0 - eta);
    N[2] = xi * eta;
    N[3] = (1.0 - xi) * eta;

    //computes the mapped point
    X[0] = N[0] * v1[0] + N[1] * v2[0] + N[2] * v3[0] + N[3] * v4[0];
    X[1] = N[0] * v1[1] + N[1] * v2[1] + N[2] * v3[1] + N[3] * v4[1];

    //computes the derivatives
    x_xi = (1.0 - eta) * (v2[0] - v1[0]) + eta * (v3[0] - v4[0]);
    x_eta = (1.0 - xi) * (v4[0] - v1[0]) + xi * (v3[0] - v2[0]);
    y_xi = (1.0 - eta) * (v2[1] - v1[1]) + eta * (v3[1] - v4[1]);
    y_eta = (1.0 - xi) * (v4[1] - v1[1]) + xi * (v3[1] - v2[1]);

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
            nodalShape(X, dX, quadraWN[1][i], quadraWN[1][j]);
            index_ij = position(i, j, q);
            Cval[index_ij][0] = (*f)(X[0], X[1]) * dX;
        }
    }

    fValSet = true;
}

double BMoment2DQuad::get_bmoment(int a1, int a2)
{
    return Bmoment[position(a1, a2, n)][0];
}

// it is expected that Fval has the function values mapped to the
// master element already multiplied by the Jacobian
// and has q * nb_Array elements
void BMoment2DQuad::setFunction(double *Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval[i][el] = Fval[i + el * q];

    fValSet = true;
}
// Fval must have at least q X nb_Array elements
void BMoment2DQuad::setFunction(double **Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval[i][el] = Fval[i][el];

    fValSet = true;
}

void BMoment2DQuad::setFunction(double (*function)(double, double))
{
    f = function;
    fDefSet = true;
}

void BMoment2DQuad::setQuadrilateral(double v1[2], double v2[2], double v3[2], double v4[2])
{
    this->v1[0] = v1[0];
    this->v1[1] = v1[1];
    this->v2[0] = v2[0];
    this->v2[1] = v2[1];
    this->v3[0] = v3[0];
    this->v3[1] = v3[1];
    this->v4[0] = v4[0];
    this->v4[1] = v4[1];
}

void BMoment2DQuad::compute_moments()
{
    if (functVal == 0 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 1 && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        if (functVal == 1)
            computeFunctionDef();

        double xi1, xi2, omega1, omega2, s1, s2, r1, r2, w1, w2;
        int i, j, a1, a2, index_ij, index_a1a2;

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                Bmoment[i][j] = 0.0;

        for (i = 0; i < q; i++)
        {
            xi1 = quadraWN[1][i];    //quadrature abcissa i
            omega1 = quadraWN[0][i]; //quadrature weight i
            s1 = 1.0 - xi1;
            r1 = xi1 / s1; //recurrence relation 1st coefficient
            for (j = 0; j < q; j++)
            {
                xi2 = quadraWN[1][j];    //quadrature abcissa j
                omega2 = quadraWN[0][j]; //quadrature weight j
                s2 = 1.0 - xi2;
                r2 = xi2 / s2; //recurrence relation 2nd coefficient
                w1 = omega1 * pow(s1, n);
                for (a1 = 0; a1 <= n; a1++)
                {
                    // here w1 equals to the weight i multiplied with the a1-th 1d bernstein polynomial evaluated at xi1
                    w2 = omega2 * pow(s2, n);
                    for (a2 = 0; a2 <= n; a2++)
                    {
                        // here w2 equals to the weight j multiplied with the a2-th 1d bernstein polynomial evaluated at xi2
                        index_ij = position(i, j, q);
                        index_a1a2 = position(a1, a2, n);
                        for (int el = 0; el < nb_Array; el ++)
                            Bmoment[index_a1a2][el] += w1 * w2 * Cval[index_ij][el];
                        w2 *= r2 * ((n - a2) / (1.0 + a2));
                    }
                    w1 *= r1 * ((n - a1) / (1.0 + a1));
                }
            }
        }
    }
}

void BMoment2DQuad::compute_moments(double (*f)(double, double))
{
    setFunction(f);
    compute_moments();
}

void BMoment2DQuad::compute_moments(double *Fval)
{
    setFunction(Fval);
    compute_moments();
}