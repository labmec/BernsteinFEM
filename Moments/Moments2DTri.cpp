#include <iostream>
#include <cmath>
#include <cstdlib>

#include "Moments.h"
#include "JacobiGaussNodes.h"

#ifdef LEN
#undef LEN
#endif
#define LEN(n, q) (MAX(n + 1, q) * MAX(n + 1, q))

BMoment2DTri::BMoment2DTri()
    : CVal(), Bmoment(), vertices(3, 2, arma::fill::none), quadraWN()
{
    std::cout << "Enter a value for the polynomial order n:";
    std::cin >> n;
    std::cout << std::endl
              << "Enter a value for the quadrature order q:";
    std::cin >> q;
    std::cout << std::endl;
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
    }

    quadraWN.set_size(q, 4);
    assignQuadra();

    lenMoments = LEN(n, q);

    Bmoment.zeros(lenMoments, 1);
    CVal.set_size(lenMoments, 1);
}

// quadrature and polynomial order constructor;
BMoment2DTri::BMoment2DTri(int q, int n)
    : CVal(LEN(n, q), 1, arma::fill::none),
      Bmoment(LEN(n, q), 1, arma::fill::zeros),
      vertices(3, 2, arma::fill::none),
      quadraWN(q, 4, arma::fill::none)
{
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
    }
    this->q = q;
    this->n = n;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);
}

// constructor setting the triangle vertices
BMoment2DTri::BMoment2DTri(int q, int n, double T[][2])
    : CVal(LEN(n, q), 1, arma::fill::none),
      Bmoment(LEN(n, q), 1, arma::fill::zeros),
      vertices(3, 2, arma::fill::none),
      quadraWN(q, 4, arma::fill::none)
{
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
    }
    this->q = q;
    this->n = n;

    assignQuadra();

    int m = MAX(n, q - 1);
    lenMoments = (m + 1) * (m + 1);

    setTriangle(T[0], T[1], T[2]);
}

BMoment2DTri::BMoment2DTri(int q, int n, int nb_Array)
    : CVal(LEN(n, q), nb_Array, arma::fill::none),
      Bmoment(LEN(n, q), nb_Array, arma::fill::zeros),
      vertices(3, 2, arma::fill::none),
      quadraWN(q, 4, arma::fill::none)
{
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
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
    : CVal(LEN(n, q), nb_Array, arma::fill::none),
      Bmoment(LEN(n, q), nb_Array, arma::fill::zeros),
      vertices(3, 2, arma::fill::none),
      quadraWN(q, 4, arma::fill::none)
{
    if (q > MAX_QUADRA_ORDER)
    {
        std::cerr << "The polynomial order is too large.\n";
        throw std::bad_alloc();
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
// function by M. Ainsworth, G. Andriamaro and O. Davydov.
void BMoment2DTri::bary2cart2d(double b1, double b2, double b3, double v1[2], double v2[2], double v3[2], double v[2])
{
    v[0] = b1 * v1[0] + b2 * v2[0] + b3 * v3[0];
    v[1] = b1 * v1[1] + b2 * v2[1] + b3 * v3[1];
}

// initialize Bmoment by the values of the function f at the quadrature points of order q
/* void BMoment2DTri::init_BmomentC_Bmom2d()
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
        b1 = quadraWN(i, 1);

        for (int j = 0; j < q; j++)
        {
            double v[2];

            index_ij = position(i, j, q - 1);

            b2 = quadraWN(j, 3) * (1 - b1);
            b3 = 1 - b1 - b2;
            //b1, b2, b3 are the barycentric coordinates using the Duffy Transformation

            bary2cart2d(b1, b2, b3, vertices.colptr(0), vertices.colptr(1), vertices.colptr(2), v); // stores b1*v1+b2*v2+b3*v3 into v;

            double functVal = (*f)(v[0], v[1]);
            Bmoment(index_ij, 0) = scalingConst * functVal;
        }
    }
}
*/

/* void BMoment2DTri::init_Bmoment2d_Cval()
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
 */
// set the function value at quadrature points, as in Fval (Fval must have at least q * q * nb_Array elements)
void BMoment2DTri::setFunction(const arma::vec &Fval)
{
    for (int i = 0; i < q * q; i++)
        for (int el = 0; el < nb_Array; el++)
            CVal(i, el) = Fval(i + el * q);

    fValSet = true;
}
// Fval must have at least q * q X nb_Array elements
void BMoment2DTri::setFunction(const arma::mat &Fval)
{
    for (int i = 0; i < q * q; i++)
        for (int el = 0; el < nb_Array; el++)
            CVal(i, el) = Fval(i, el);

    fValSet = true;
}

// set the function definition for computations
void BMoment2DTri::setFunction(std::function<double(double, double)> function)
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

void BMoment2DTri::setTriangle(const arma::mat &vertices)
{
    this->vertices = vertices;
}

void BMoment2DTri::computeFunctionDef()
{
    arma::mat points(getIntegrationPoints());

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            int index_ij = position_q(i, j, q);
            CVal(index_ij, 0) = f(points(i, 0), points(j, 1));
        }
    }
    fValSet = true;
}

// returns the vectors with the integration points (x, y) over the object's element, following the moments organization
// Assuming: points = getIntegrationPoints(); then
// points(i, 0) == x i-th coordinate
// points(i, 1) == y i-th coordinate
arma::mat BMoment2DTri::getIntegrationPoints()
{
    arma::mat points(q * q * nb_Array, 2); // vector with q * q * nb_Array elements
    double b1, b2, b3;
    double v1[2] = {vertices(0, 0), vertices(0, 1)};
    double v2[2] = {vertices(1, 0), vertices(1, 1)};
    double v3[2] = {vertices(2, 0), vertices(2, 1)};

    for (int i = 0; i < q; i++)
    {
        b1 = quadraWN(i, 1);
        for (int j = 0; j < q; j++)
        {
            double v[2];
            b2 = quadraWN(j, 3) * (1 - b1);
            b3 = 1 - b1 - b2;
            bary2cart2d(b1, b2, b3, v1, v2, v3, v); // stores b1*v1+b2*v2+b3*v3 into v;
            int index_ij = position_q(i, j, q);
            points(index_ij, 0) = v[0];
            points(index_ij, 1) = v[1];
        }
    }

    return points;
}

//compute the b-moments
void BMoment2DTri::compute_moments()
{
    if (functVal == 1 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 0 && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        int m = MAX(n, q - 1); // m will be used for indexing

        arma::mat Bmoment_inter((m + 1) * (m + 1), nb_Array, arma::fill::zeros);
        Bmoment.zeros();

        // compute the function definition into the function values vector
        if (functVal == 0)
            computeFunctionDef();

        double scalingConst = 2 * Area2d(vertices); // NEED this scaling constant, as result depends on Area2d(v1,v2,v3)
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
                    int index_a1j = position_q(a1, j, m);

                    int index_ij = position_q(i, j, q);

                    for (int ell = 0; ell < nb_Array; ell++)
                    {
                        Bmoment_inter(index_a1j, ell) += scalingConst * B * CVal(index_ij, ell);
                    }
                }
                B = B * r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=1)
        for (int i = 0; i < q; i++)
        {
            double xi = quadraWN(i, 3);
            double wgt = quadraWN(i, 2);

            double s = 1 - xi;
            double r = xi / s;

            for (int a1 = 0; a1 <= n; a1++)
            {
                double B = wgt * pow(s, n - a1);
                for (int a2 = 0; a2 <= n - a1; a2++)
                {
                    int index_a1a2 = position(a1, a2, n);
                    int index_a1i = position_q(a1, i, m);

                    for (int ell = 0; ell < nb_Array; ell++)
                    {
                        Bmoment(index_a1a2, ell) += B * Bmoment_inter(index_a1i, ell);
                    }
                    B = B * r * (n - a1 - a2) / (1 + a2);
                }
            }
        }
    }
}

void BMoment2DTri::compute_moments(std::function<double(double, double)> f)
{
    setFunction(f);
    compute_moments();
}

void BMoment2DTri::compute_moments(const arma::vec &Fval)
{
    setFunction(Fval);
    compute_moments();
}

/* modified from 'bbfem.cpp' to fit into this implementation */
// Deprecated
void BMoment2DTri::transform_BmomentC_Stiff2d(BMoment2DTri *Bmomentab, const arma::mat &normalMat)
{
    int m, mm;

    //m = MAX (2 * n - 2, q-1);
    m = MAX(n, q - 1); // n from this object is already expected to be (2 * n - 2) from the StiffM2DTri object

    mm = (m + 1) * (m + 1); // with no storing algorithm, required Bmoment size is determined by MAX( Bmoment order, number of quadrature points )

    for (int mu = 0; mu < mm; mu++)
    {
        arma::mat Mat(2, 2, arma::fill::none);

        Mat.at(0, 0) = Bmoment.at(mu, 0); // Mat is used to store Bmoment[mu] entries
        Mat.at(0, 1) = Mat.at(1, 0) = Bmoment.at(mu, 1);
        Mat.at(1, 1) = Bmoment.at(mu, 2);

        double matVectMult[2];

        matVectMult[0] = Mat.at(0, 0) * normalMat.at(0, 0) + Mat.at(0, 1) * normalMat.at(0, 1);
        matVectMult[1] = Mat.at(1, 0) * normalMat.at(0, 0) + Mat.at(1, 1) * normalMat.at(0, 1);

        Bmomentab->Bmoment.at(mu, 0) = normalMat.at(0, 0) * matVectMult[0] + normalMat.at(0, 1) * matVectMult[1]; //alfa=[1,0,0], beta=[1,0,0]
        Bmomentab->Bmoment.at(mu, 1) = normalMat.at(1, 0) * matVectMult[0] + normalMat.at(1, 1) * matVectMult[1]; //alfa=[0,1,0], beta=[1,0,0]
        Bmomentab->Bmoment.at(mu, 2) = normalMat.at(2, 0) * matVectMult[0] + normalMat.at(2, 1) * matVectMult[1]; //alfa=[0,0,1], beta=[1,0,0]

        matVectMult[0] = Mat.at(0, 0) * normalMat.at(1, 0) + Mat.at(0, 1) * normalMat.at(1, 1);
        matVectMult[1] = Mat.at(1, 0) * normalMat.at(1, 0) + Mat.at(1, 1) * normalMat.at(1, 1);

        //Bmomentab->Bmoment.at(mu, 1) = normalMat.at(0, 0)*matVectMult[0] + normalMat.at(0, 1)*matVectMult[1]; //alfa=[1,0,0], beta=[0,1,0] // redundancy
        Bmomentab->Bmoment.at(mu, 3) = normalMat.at(1, 0) * matVectMult[0] + normalMat.at(1, 1) * matVectMult[1]; //alfa=[0,1,0], beta=[0,1,0]
        Bmomentab->Bmoment.at(mu, 4) = normalMat.at(2, 0) * matVectMult[0] + normalMat.at(2, 1) * matVectMult[1]; //alfa=[0,0,1], beta=[0,1,0]

        matVectMult[0] = Mat.at(0, 0) * normalMat.at(2, 0) + Mat.at(0, 1) * normalMat.at(2, 1);
        matVectMult[1] = Mat.at(1, 0) * normalMat.at(2, 0) + Mat.at(1, 1) * normalMat.at(2, 1);

        //Bmomentab->Bmoment.at(mu, 2) = normalMat.at(0, 0)*matVectMult[0] + normalMat.at(0, 1)*matVectMult[1]; //alfa=[1,0,0], beta=[0,0,1] // redundancy
        //Bmomentab->Bmoment.at(mu. 4) = normalMat.at(1, 0)*matVectMult[0] + normalMat.at(1, 1)*matVectMult[1]; //alfa=[0,1,0], beta=[0,0,1]
        Bmomentab->Bmoment.at(mu, 5) = normalMat.at(2, 0) * matVectMult[0] + normalMat.at(2, 1) * matVectMult[1]; //alfa=[0,0,1], beta=[0,0,1]
    }
}