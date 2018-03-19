#include <iostream>
using std::cout;
using std::cin;
using std::endl;

#include "Moments.h"
#include "JacobiGaussNodes.h"

// default constructor
BMoment1D::BMoment1D ()
{
    cout << "Enter a value for the polynomial order n:";
    cin >> n;
    cout << endl << "Enter a value for the quadrature order q:";
    cin >> q;
    cout << endl;
    if (q>80)
    {
	    std::cerr<<"The polynomial order is too large.\n";
	    exit (EXIT_FAILURE);
    }

    quadraWN = create_quadraWN();
    assignQuadra();

    int m = MAX(n, q-1);
    lenMoments = (m+1);

    Bmoment = create_Bmoment();
}

// quadrature and polynomial order constructor
BMoment1D::BMoment1D (int q, int n)
{
    if (q>80)
    {
	    std::cerr<<"The polynomial order is too large.\n";
	    exit (EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;

    quadraWN = create_quadraWN();
    assignQuadra();

    int m = MAX(n, q-1);
    lenMoments = (m+1);

    Bmoment = create_Bmoment();
}

BMoment1D::~BMoment1D ()
{
    if(fValSet)
        delete Cval;
    delete Bmoment;
    delete quadraWN[0];
    delete quadraWN;
}


// alloc Bmoment vector
double* BMoment1D::create_Bmoment ()
{
    return new double[lenMoments];
}

// free Bmoment vector memory
void BMoment1D::delete_Bmoment (double *Bmoment) { delete Bmoment; }

// alloc memory for the quadrature points and weights
double** BMoment1D::create_quadraWN ()
{
    double *aux = new double [2*(q+1)];
    double **quadraWN = new double*[2];
    quadraWN[0] = aux;
    quadraWN[1] = aux + q + 1;
}

void BMoment1D::assignQuadra ()
{        
    double *x = quadraWN[1];
    double *w = quadraWN[0];

    for (int k = 0; k < q; k++)
    {
        x[k] = ( 1.0 + legendre[1][q-2][k] ) * 0.5;
        w[k] = legendre[0][q-2][k] * 0.5;
    }
}

// set the function values for computation
void BMoment1D::setFunctionValue (double *Fval)
{
    Cval = new double[q];
    for (int i = 0; i < q; i++)
        Cval[i] = Fval[i];
    fValSet = true;
}

// set the function definition for computation
void BMoment1D::setFunction (double (*function) (double))
{
    f = function;
    fDefSet = true;
}

// compute the B-moments using the values already assigned in the object
void BMoment1D::compute_moments ()
{
    if(functVal == 0 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 1 && !fValSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else {
        int i, alpha;
        double xi, omega, s, r, w;

        for (i = 0; i <= n; Bmoment[i] = 0, i++); // sets Bmoment to 0 in all entries

        for (i = 0; i < q; i++) {
            xi = quadraWN[1][i];
            omega = quadraWN[0][i];
            s = 1 - xi;
            r = xi / s;
            w = omega * pow(s, n);
            for (alpha = 0; alpha <= n; alpha++) {
                // here w equals the Bernstein polynomial of order n,
                // with index alpha, evaluated at the i-th integration node
                // times the i-th integration weight.
                Bmoment[alpha] += w * Cval[i];
                w *= r * ((n - alpha) / (1. + alpha)); // treats the recurrence relation
            }
        }
    }
}

    // compute the b-moments for the specified f function
    void BMoment1D::compute_moments (double (*f) (double))
{
    setFunction(f);
    compute_moments();
}

// compute the b-moments for the Fval function values
void BMoment1D::compute_moments (double *Fval)
{
    setFunctionValue(Fval);
    compute_moments();
}