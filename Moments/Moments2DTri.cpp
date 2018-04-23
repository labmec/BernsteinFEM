#include <iostream>
#include <cmath>
#include <cstdlib>

#include "Moments.h"
#include "JacobiGaussNodes.h"

BMoment2DTri::BMoment2DTri()
    : Bmoment(), CVal(), quadraWN(), vertices(3, 2, arma::fill::none)
{
    std::cout << "Enter a value for the polynomial order n:";
    std::cin >> n;
    std::cout << std::endl
              << "Enter a value for the quadrature order q:";
    std::cin >> q;
    std::cout << std::endl;
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }

    quadraWN.set_size(q, 4);
    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);

    Bmoment.zeros(lenMoments, 1);
    CVal.set_size(q, 1);
}

// quadrature and polynomial order constructor;
BMoment2DTri::BMoment2DTri(int q, int n)
    : Bmoment(MAX(n + 1, q) * MAX(n + 1, q), 1, arma::fill::zeros),
      CVal(q * q, 1, arma::fill::none),
      quadraWN(q, 4, arma::fill::none),
      vertices(3, 2, arma::fill::none)
{
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);
}

// constructor setting the triangle vertices
BMoment2DTri::BMoment2DTri(int q, int n, double T[][2])
    : Bmoment(MAX(n + 1, q) * MAX(n + 1, q), 1, arma::fill::zeros),
      CVal(q * q, 1, arma::fill::none),
      quadraWN(q, 4, arma::fill::none),
      vertices(3, 2, arma::fill::none)
{
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);

    setTriangle(T[0], T[1], T[2]);
}

BMoment2DTri::BMoment2DTri(int q, int n, int nb_Array)
    : Bmoment(MAX(n + 1, q) * MAX(n + 1, q), nb_Array, arma::fill::zeros),
      CVal(q * q, nb_Array, arma::fill::none),
      quadraWN(q, 4, arma::fill::none),
      vertices(3, 2, arma::fill::none)
{
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;
    this->nb_Array = nb_Array;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);
}

// constructor setting the triangle vertices
BMoment2DTri::BMoment2DTri(int q, int n, int nb_Array, double T[][2])
    : Bmoment(MAX(n + 1, q) * MAX(n + 1, q), nb_Array, arma::fill::zeros),
      CVal(q * q, nb_Array, arma::fill::none),
      quadraWN(q, 4, arma::fill::none),
      vertices(3, 2, arma::fill::none)
{
    if (q > 80)
    {
        std::cerr << "The polynomial order is too large.\n";
        exit(EXIT_FAILURE);
    }
    this->q = q;
    this->n = n;
    this->nb_Array = nb_Array;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);

    setTriangle(T[0], T[1], T[2]);
}

BMoment2DTri::~BMoment2DTri()
{
}

void BMoment2DTri::assignQuadra()
{
    int k;

    double *x = quadraWN.colptr(1); // x is pointer which points to the address l_x, thus the effective MODIFICATION of l_x entries
    double *w = quadraWN.colptr(0);

    for (k = 0; k < q; k++)
    {
        x[k] = (1.0 + jacobi_xi(q, k)) * 0.5;
        w[k] = jacobi_w(q, k) * 0.25;
    }

    x = quadraWN.colptr(3);
    w = quadraWN.colptr(2);

    for (k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

//convert barycentric coordinates (b1,b2,b3) of a point w.r.t. vertices v1,v2,v3 into Cartesian coordinates v
void BMoment2DTri::bary2cart2d(double b1, double b2, double b3, double v1[2], double v2[2], double v3[2], double v[2])
{
    v[0] = b1 * v1[0] + b2 * v2[0] + b3 * v3[0];
    v[1] = b1 * v1[1] + b2 * v2[1] + b3 * v3[1];
}

// initialize Bmoment by the values of the function f at the quadrature points of order q
void BMoment2DTri::init_BmomentC_Bmom2d()
{
    if (nb_Array > 1)
    {
        std::cerr << "Function definition may be used only scalar-valued" << std::endl;
        nb_Array = 1;
    }
    double b1, b2, b3;
    int index_ij = 0;

    double scalingConst = 2 * Area2d(vertices); // NEED this scaling constant, as result depends on Area2d(v1,v2,v3)

    for (int i = 0; i < q; i++)
    {
        b1 = quadraWN(1, i);

        for (int j = 0; j < q; j++)
        {
            double v[2];

            index_ij = position(i, j, q - 1);

            b2 = quadraWN(3, j) * (1 - b1);
            b3 = 1 - b1 - b2;
            //b1, b2, b3 are the barycentric coordinates using the Duffy Transformation

            bary2cart2d(b1, b2, b3, vertices.colptr(0), vertices.colptr(1), vertices.colptr(2), v); // stores b1*v1+b2*v2+b3*v3 into v;

            double functVal = (*f)(v[0], v[1]);
            Bmoment(index_ij, 0) = scalingConst * functVal;
        }
    }
}

void BMoment2DTri::init_Bmoment2d_Cval()
{
    double scalingConst = 2 * Area2d(vertices); // Jacobian of Duffy transformation

    int index_ij = 0;

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int ell = 0; ell < nb_Array; ell++)
                Bmoment(index_ij, ell) = scalingConst * CVal(index_ij, ell);
            index_ij++;
        }
    }
}

//get the bmoment value of the Bernstein polynomial with indexes a1 and a2 (a3 = n - a2 - a1) on the specified dimension
double BMoment2DTri::get_bmoment(int a1, int a2, int dim)
{
    int i = position(a1, a2, n);
    return Bmoment(i, dim - 1);
}

// set the function value at quadrature points, as in Fval (Fval must have at least q * nb_Array elements)
void BMoment2DTri::setFunction(arma::vec Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            CVal(i, el) = Fval[i + el * q];

    fValSet = true;
}
// Fval must have at least q X nb_Array elements
void BMoment2DTri::setFunction(arma::mat Fval)
{
    for (int i = 0; i < q; i++)
        for (int el = 0; el < nb_Array; el++)
            CVal(i, el) = Fval(i, el);

    fValSet = true;
}

// set the function definition for computations
void BMoment2DTri::setFunction(double (*function)(double, double))
{
    f = function;
    fDefSet = true;
}

// set the element triangle vertices
void BMoment2DTri::setTriangle(double v1[2], double v2[2], double v3[2])
{
    vertices(0, 0) = v1[0];
    vertices(0, 1) = v1[1];
    vertices(1, 0) = v2[0];
    vertices(1, 1) = v2[1];
    vertices(2, 0) = v3[0];
    vertices(2, 1) = v3[1];
}

//compute the b-moments
void BMoment2DTri::compute_moments()
{
    if (functVal == 0 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 1 && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        int m = MAX(n, q - 1); // m will be used for indexing

        arma::mat BmomentInter((m + 1) * (m + 1), 1, arma::fill::zeros);

        //initialize Bmoment with function values
        if (functVal == 1)
            if (fValSet)
                init_Bmoment2d_Cval();
            else
                std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
        else if (fDefSet)
            init_BmomentC_Bmom2d();
        else
            std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";

        //  initialize BmomentInter, at this point, array Bmoment has been initialized with function values
        for (int i = 0; i <= m; i++)
        {
            for (int j = 0; j <= m; j++)
            {
                int index_ij = position(i, j, m);

                for (int ell = 0; ell < nb_Array; ell++)
                    BmomentInter(index_ij, ell) = 0;
            }
        }

        double xi, wgt, s, r, B;

        arma::vec fact(nb_Array, arma::fill::none);

        // convert first index (l=2)
        for (int i = 0; i < q; i++)
        {
            xi = quadraWN(1, i);
            wgt = quadraWN(0, i);

            s = 1 - xi;
            r = xi / (1 - xi);

            B = wgt * pow(s, n);
            for (int mu1 = 0; mu1 <= n; mu1++)
            {
                for (int j = 0; j < q; j++)
                {
                    int index_mu1j = position(mu1, j, m);

                    int index_ij = position(i, j, q - 1); // init_BmomentC uses this indexing

                    for (int ell = 0; ell < nb_Array; ell++)
                    {
                        fact(ell) = B * Bmoment(index_ij, ell);

                        BmomentInter(index_mu1j, ell) += fact(ell);
                    }
                }
                B *= r * (n - mu1) / (1 + mu1);
            }
        }

        // reset the Bmoment array
        for (int mu1 = 0; mu1 <= n; mu1++)
        {
            for (int mu2 = 0; mu2 <= n - mu1; mu2++)
            {
                int index_mu1mu2 = position(mu1, mu2, m);

                for (int ell = 0; ell < nb_Array; ell++)
                    Bmoment(index_mu1mu2, ell) = 0;
            }
        }

        // convert second index (l=1)
        for (int j = 0; j < q; j++)
        {
            xi = quadraWN(3, j);
            wgt = quadraWN(2, j);

            s = 1 - xi;
            r = xi / (1 - xi);

            for (int mu1 = 0; mu1 <= n; mu1++)
            {
                B = wgt * pow(s, n - mu1);
                for (int mu2 = 0; mu2 <= n - mu1; mu2++)
                {
                    int index_mu1mu2 = position(mu1, mu2, m);
                    int index_mu1j = position(mu1, j, m);

                    for (int ell = 0; ell < nb_Array; ell++)
                    {
                        fact(ell) = B * BmomentInter(index_mu1j, ell);
                        Bmoment(index_mu1mu2, ell) += fact(ell);
                    }
                    B *= r * (n - mu1 - mu2) / (1 + mu2);
                }
            }
        }
    }
}

void BMoment2DTri::compute_moments(double (*f)(double, double))
{
    setFunction(f);
    compute_moments();
}

void BMoment2DTri::compute_moments(arma::vec Fval)
{
    setFunction(Fval);
    compute_moments();
}

/* modified from 'bbfem.cpp' to fit into this implementation */
void BMoment2DTri::transform_BmomentC_Stiff2d(BMoment2DTri *Bmomentab, arma::mat normalMat)
{
    int m, mm;

    //m = MAX (2 * n - 2, q-1);
    m = MAX(n, q - 1); // n from this object is already expected to be (2 * n - 2) from the StiffM2DTri object

    mm = (m + 1) * (m + 1); // with no storing algorithm, required Bmoment size is determined by MAX( Bmoment order, number of quadrature points )

    for (int mu = 0; mu < mm; mu++)
    {
        arma::mat Mat(2, 2, arma::fill::none);

        Mat(0, 0) = Bmoment(mu, 0); // Mat is used to store Bmoment[mu] entries
        Mat(0, 1) = Mat(1, 0) = Bmoment(mu, 1);
        Mat(1, 1) = Bmoment(mu, 2);

        double matVectMult[2];

        matVectMult[0] = Mat(0, 0) * normalMat(0, 0) + Mat(0, 1) * normalMat(0, 1);
        matVectMult[1] = Mat(1, 0) * normalMat(0, 0) + Mat(1, 1) * normalMat(0, 1);

        Bmomentab->Bmoment(mu, 0) = normalMat(0, 0) * matVectMult[0] + normalMat(0, 1) * matVectMult[1]; //alfa=[1,0,0], beta=[1,0,0]
        Bmomentab->Bmoment(mu, 1) = normalMat(1, 0) * matVectMult[0] + normalMat(1, 1) * matVectMult[1]; //alfa=[0,1,0], beta=[1,0,0]
        Bmomentab->Bmoment(mu, 2) = normalMat(2, 0) * matVectMult[0] + normalMat(2, 1) * matVectMult[1]; //alfa=[0,0,1], beta=[1,0,0]

        matVectMult[0] = Mat(0, 0) * normalMat(1, 0) + Mat(0, 1) * normalMat(1, 1);
        matVectMult[1] = Mat(1, 0) * normalMat(1, 0) + Mat(1, 1) * normalMat(1, 1);

        //Bmomentab->Bmoment(mu, 1) = normalMat(0, 0)*matVectMult[0] + normalMat(0, 1)*matVectMult[1]; //alfa=[1,0,0], beta=[0,1,0] // redundancy
        Bmomentab->Bmoment(mu, 3) = normalMat(1, 0) * matVectMult[0] + normalMat(1, 1) * matVectMult[1]; //alfa=[0,1,0], beta=[0,1,0]
        Bmomentab->Bmoment(mu, 4) = normalMat(2, 0) * matVectMult[0] + normalMat(2, 1) * matVectMult[1]; //alfa=[0,0,1], beta=[0,1,0]

        matVectMult[0] = Mat(0, 0) * normalMat(2, 0) + Mat(0, 1) * normalMat(2, 1);
        matVectMult[1] = Mat(1, 0) * normalMat(2, 0) + Mat(1, 1) * normalMat(2, 1);

        //Bmomentab->Bmoment(mu, 2) = normalMat(0, 0)*matVectMult[0] + normalMat(0, 1)*matVectMult[1]; //alfa=[1,0,0], beta=[0,0,1] // redundancy
        //Bmomentab->Bmoment(mu. 4) = normalMat(1, 0)*matVectMult[0] + normalMat(1, 1)*matVectMult[1]; //alfa=[0,1,0], beta=[0,0,1]
        Bmomentab->Bmoment(mu, 5) = normalMat(2, 0) * matVectMult[0] + normalMat(2, 1) * matVectMult[1]; //alfa=[0,0,1], beta=[0,0,1]
    }
}