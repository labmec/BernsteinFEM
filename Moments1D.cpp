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
    lenMoments = (m + 1);

    Bmoment = create_Bmoment();
    Cval = create_Cval();
}

// quadrature and polynomial order constructor
BMoment1D::BMoment1D(int q, int n)
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
    lenMoments = (m + 1);

    Bmoment = create_Bmoment();
    Cval = create_Cval();
}

BMoment1D::~BMoment1D()
{
    delete Cval[0];
    delete Cval;
    delete Bmoment;
    delete quadraWN[0];
    delete quadraWN;
}

// alloc Bmoment vector
double **BMoment1D::create_Bmoment()
{
    double *aux = new double[lenMoments * nb_Array];
    double **bmoment = new double *[lenMoments];

    for (int i = 0; i < lenMoments; aux += nb_Array, i++)
        bmoment[i] = aux;

    return bmoment;
}

double **BMoment1D::create_Cval()
{
    int aux = lenMoments;
    lenMoments = q;
    double **p = create_Bmoment();
    lenMoments = aux;
    return p;
}

// free Bmoment vector memory
void BMoment1D::delete_Bmoment(double **Bmoment)
{
    delete Bmoment[0];
    delete Bmoment;
}

// alloc memory for the quadrature points and weights
double **BMoment1D::create_quadraWN()
{
    double *aux = new double[2 * (q + 1)];
    double **quadraWN = new double *[2];
    quadraWN[0] = aux;
    quadraWN[1] = aux + q + 1;
}

void BMoment1D::assignQuadra()
{
    double *x = quadraWN[1];
    double *w = quadraWN[0];

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

void BMoment1D::loadFunctionDef()
{
    if (nb_Array > 1)
    {
        std::cerr << "Function definition may be used only scalar-valued" << endl;
        nb_Array = 1;
    }
    double *aux = new double[q * nb_Array];
    Cval = new double *[q];

    for (int i = 0; i < q; aux += nb_Array, i++)
    {
        Cval[i] = aux;
        for (int el = 0; el < nb_Array; el++)
            Cval[i][el] = (*f)(quadraWN[0][i]);
    }

    fValSet = true;
}

// set the function values for computation (Fval must have at least q * nb_Array elements)
void BMoment1D::setFunction(double *Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval[i][el] = Fval[i + el * q];

    fValSet = true;
}
// Fval must have at least q X nb_Array elements
void BMoment1D::setFunction(double **Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            Cval[i][el] = Fval[i][el];

    fValSet = true;
}

// set the function definition for computation
void BMoment1D::setFunction(double (*function)(double))
{
    f = function;
    fDefSet = true;
}

void BMoment1D::setInterval(double a, double b)
{
    this->a = a;
    this->b = b;
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
        int i, alpha;
        double xi, omega, s, r, w;

        for (i = 0; i <= n; Bmoment[i] = 0, i++)
            ; // sets Bmoment to 0 in all entries

        for (i = 0; i < q; i++)
        {
            xi = quadraWN[1][i];
            omega = quadraWN[0][i];
            s = 1 - xi;
            r = xi / s;
            w = omega * pow(s, n);
            for (alpha = 0; alpha <= n; alpha++)
            {
                // here w equals the Bernstein polynomial of order n,
                // with index alpha, evaluated at the i-th integration node
                // times the i-th integration weight.
                for (int el = 0; el < nb_Array; el++)
                    Bmoment[alpha][el] += w * Cval[i][el];
                w *= r * ((n - alpha) / (1. + alpha)); // treats the recurrence relation
            }
        }
    }
}

// compute the b-moments for the specified f function
void BMoment1D::compute_moments(double (*f)(double))
{
    setFunction(f);
    compute_moments();
}

// compute the b-moments for the Fval function values
void BMoment1D::compute_moments(double *Fval)
{
    setFunction(Fval);
    compute_moments();
}