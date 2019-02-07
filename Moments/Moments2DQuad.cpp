#include "Moments.h"
#include "JacobiGaussNodes.h"

// helps indexing quadrature points vectors
uint position_q(uint i, uint j, int q) { return i * q + j; }

BMoment2DQuad::BMoment2DQuad(uint q, uint n, const Element<Element_t::QuadrilateralEl> &element, uint nb_Array)
    : BMoment2DQuad(q, n, n, element, nb_Array) {}

BMoment2DQuad::BMoment2DQuad(uint q, uint n, uint m, const Element<Element_t::QuadrilateralEl> &element, uint nb_Array)
    : BMoment(q, n, element, nb_Array),
      BMomentInter((MAX(n + 1, q)) * (MAX(n + 1, q)), nb_Array, arma::fill::zeros)
{
    this->m = m;
    int max = MAX(n, m);
    lenMoments = (max + 1) * (max + 1);
    lenCval = q * q;
    Bmoment.set_size(lenMoments, nb_Array);
    quadraWN.set_size(q, 2);
    Cval.set_size(lenCval, nb_Array);
    assignQuadra();
}

// copy constructor
BMoment2DQuad::BMoment2DQuad(const BMoment2DQuad &cp)
    : BMoment(cp.q, cp.n, cp.element, cp.nb_Array),
      BMomentInter(cp.BMomentInter)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
    quadraWN = cp.quadraWN;
}

// copy assignment operator
BMoment2DQuad &BMoment2DQuad::operator=(const BMoment2DQuad &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        nb_Array = cp.nb_Array;
        lenMoments = cp.lenMoments;
        lenCval = cp.lenCval;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
        quadraWN = cp.quadraWN;
    }
    return *this;
}

inline
void BMoment2DQuad::assignQuadra()
{
    double *x = quadraWN.colptr(1);
    double *w = quadraWN.colptr(0);

    for (uint k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment2DQuad::loadFunctionDef()
{

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            arma::mat jac(2, 2, arma::fill::none);
            arma::mat X = element.mapToElement({quadraWN.at(i, 1), quadraWN.at(j, 1)}, jac);
            double dX = det(jac);
            uint index_ij = position_q(i, j, q);
            for (uint el = 0; el < nb_Array; el++)
                Cval.at(index_ij, 0) = f(X[0], X[1]) * dX;
        }
    }

    fValSet = true;
}

inline
arma::mat BMoment2DQuad::getIntegrationPoints()
{
    arma::mat points(lenCval, 2, arma::fill::none);

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            arma::mat X = element.mapToElement({quadraWN.at(i, 1), quadraWN.at(j, 1)});
            uint index_ij = position_q(i, j, q);
            points.at(index_ij, 0) = X[0];
            points.at(index_ij, 1) = X[1];
        }
    }

    return points;
}

arma::mat &BMoment2DQuad::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        if (!functVal)
            loadFunctionDef();

        uint max_nm = MAX(n, m);
        uint max_nq = MAX(max_nm, q - 1);
        Bmoment.zeros();
        BMomentInter.zeros();

        // convert first index (l=2)
        for (uint i = 0; i < q; i++)
        {
            double xi = quadraWN(i, 1);
            double wgt = quadraWN(i, 0);

            double s = 1 - xi;
            double r = xi / s;

            double B = wgt * pow(s, n);
            for (uint a1 = 0; a1 <= n; a1++)
            {
                for (uint j = 0; j < q; j++)
                {
                    uint index_a1j = position_q(a1, j, max_nq);

                    uint index_ij = position_q(i, j, q);

                    for (uint ell = 0; ell < nb_Array; ell++)
                        BMomentInter(index_a1j, ell) += B * Cval(index_ij, ell);
                }
                B = B * r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=1)
        for (uint i = 0; i < q; i++)
        {
            double xi = quadraWN(i, 1);
            double wgt = quadraWN(i, 0);

            double s = 1 - xi;
            double r = xi / s;

            double B = wgt * pow(s, m);
            for (uint a2 = 0; a2 <= m; a2++)
            {
                for (uint a1 = 0; a1 <= n; a1++)
                {
                    uint index_a1i = position_q(a1, i, max_nq);
                    uint index_a1a2 = element.position({a1, a2}, max_nm);

                    for (uint ell = 0; ell < nb_Array; ell++)
                        Bmoment(index_a1a2, ell) += B * BMomentInter(index_a1i, ell);
                }
                B = B * r * (m - a2) / (1 + a2);
            }
        }
    }
    return Bmoment;
}
