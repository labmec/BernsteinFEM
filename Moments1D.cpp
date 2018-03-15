#include <iostream>
#include "Moments.h"
#include "JacobiGaussNodes.h"

class BMoment1D
{
    int lenMoments; // length of the Bmoment vector
    double *Bmoment; // vector where the b-moments are stored
    double *Cval; // vector where the function values are stored
    bool fValSet = false; // is true if the function value is set
    bool fDefSet = false; // is true if the function definition is set

    // alloc Bmoment vector
    double *create_Bmoment()
    {
        return new double[lenMoments];
    }

    // free Bmoment vector memory
    void delete_Bmoment(double *Bmoment) { delete Bmoment; }

    // alloc memory for the quadrature points and weights
    double** create_quadraWN ()
    {
        double *aux = new double [2*(q+1)];
        double **quadraWN = new double*[2];
        quadraWN[0] = aux;
        quadraWN[1] = aux + q + 1;
    }

    void assignQuadra()
    {        
        double *x = quadraWN[1];
        double *w = quadraWN[0];

        for (int k = 0; k < q; k++)
        {
            x[k] = ( 1.0 + legendre[1][q-2][k] ) * 0.5;
            w[k] = legendre[0][q-2][k] * 0.5;
        }
    }

protected:

    int q; // number of quadrature points in one dimension 
    int n; // Bernstein polynomial order
    int functVal = 1; // determines if will use function value or function definition (use function dvalue by default)
    int nb_Array = 1; // dimension of function Image (function is scalar valued by default)
    double a, b; // interval [a, b] variables
    double **quadraWN; // quadrature points and weights

    // function definition for the computation of the b-moments
    double (*f) (double) = 0x0;

public:
    // default constructor
    BMoment1D()
    {
        std::cout << "Enter a value for the polynomial order n:";
        std::cin >> n;
        std::cout << std::endl << "Enter a value for the quadrature order q:";
        std::cin >> q;
        std::cout << std::endl;
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
    BMoment1D(int q, int n)
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

    ~BMoment1D()
    {
        if(fValSet)
            delete Cval;
        delete Bmoment;
        delete quadraWN[0];
        delete quadraWN;
    }

    // set the function values for computation
    void setFunctionValue (double *Fval)
    {
        Cval = new double[q];
        for (int i = 0; i < q; i++)
            Cval[i] = Fval[i];
        fValSet = true;
    }

    // set the function definition for computation
    void setFunction (double (*function) (double))
    {
        f = function;
        fDefSet = true;
    }

    // compute the B-moments using the values already assigned in the object
    void compute_moments ()
    {
         if(!fValSet && ! fDefSet)
            std::cerr << "missing function values or definition for computation of the moments in \'compute_moments()\'\n";
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
    void compute_moments (double (*f) (double))
    {
        setFunction(f);
        compute_moments();
    }

    // compute the b-moments for the Fval function values
    void compute_moments (double *Fval)
    {
        setFunctionValue(Fval);
        compute_moments();
    }
};