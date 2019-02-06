#include "Moments.h"
#include "JacobiGaussNodes.h"

uint position_q(uint i, uint j, uint k, int q) { return i * q * q + j * q + k; }

BMoment3DCube::BMoment3DCube(int q, int n, const Element<Element_t::CubeEl> &el, int nb_Array)
    : BMoment3DCube(q, n, n, n, el, nb_Array) {}

BMoment3DCube::BMoment3DCube(int q, int n, int m, int p, const Element<Element_t::CubeEl> &el, int nb_Array)
    : BMoment(q, n, el, nb_Array), 
    BMomentInter()
{
    int max = MAX(MAX(n, m), p);
    max++;
    lenMoments = max * max * max; // (max + 1) ^ 3
    lenCval = q * q * q;
    Bmoment.set_size(lenMoments, nb_Array);
    BMomentInter.set_size(lenMoments, nb_Array);
    quadraWN.set_size(q, 2);
    Cval.set_size(lenCval, nb_Array);
    assignQuadra();
}

// copy constructor
BMoment3DCube::BMoment3DCube(const BMoment3DCube &cp)
    : BMoment(cp.q, cp.n, cp.element, cp.nb_Array), BMomentInter(cp.BMomentInter)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
    quadraWN = cp.quadraWN;
}

BMoment3DCube &BMoment3DCube::operator=(const BMoment3DCube &cp)
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
void BMoment3DCube::assignQuadra()
{
    double *x = quadraWN.colptr(1);
    double *w = quadraWN.colptr(0);

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment3DCube::loadFunctionDef()
{
    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; i++)
        {
            for (int k = 0; k < q; k++)
            {
                arma::mat jac(3, 3, arma::fill::none);
                arma::mat X = element.mapToElement({quadraWN(i, 1), quadraWN(j, 1), quadraWN(k, 1)}, jac);
                double dX = det(jac);
                uint index_ijk = position_q(i, j, k, q);
                for (int el = 0; el < nb_Array; el++)
                    Cval.at(index_ijk, el) = f(X[0], X[1], X[2]) * dX;
            }
        }
    }
}

inline
arma::mat BMoment3DCube::getIntegrationPoints()
{
    arma::mat points(lenCval, 3, arma::fill::none);

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < q; k++)
            {
                arma::mat X = element.mapToElement({quadraWN(i, 1), quadraWN(j, 1), quadraWN(k, 1)});
                uint index_ijk = position_q(i, j, k, q);
                points.at(index_ijk, 0) = X[0];
                points.at(index_ijk, 1) = X[1];
                points.at(index_ijk, 2) = X[2];
            }
        }
    }

    return points;
}

arma::mat &BMoment3DCube::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        if (!functVal)
            loadFunctionDef();

        int max_nm = MAX(MAX(n, m), p);
        int max_nq = MAX(max_nm, q);
        Bmoment.zeros();
        BMomentInter.zeros();

        // convert first index (l=3)
        for (uint i = 0; i < q; i++)
        {
            double xi = quadraWN.at(i, 1);
            double wgt = quadraWN.at(i, 0);

            double s = 1 - xi;
            double r = xi / s;
            double B = wgt * pow(s, n);
            for (uint a1 = 0; a1 <= n; a1++)
            {
                for (uint j = 0; j < q; j++)
                {
                    for (uint k = 0; j < q; k++)
                    {
                        uint index_a1jk = position_q(a1, j, k, max_nq);
                        uint index_ijk = position_q(i, j, k, q);

                        for (uint el = 0; el < nb_Array; el++)
                            Bmoment.at(index_a1jk, el) += B * Cval.at(index_ijk, el);
                    }
                }
                B *= r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=2)
        for (uint j = 0; j < q; j++)
        {
            double xi = quadraWN.at(j, 1);
            double wgt = quadraWN.at(j, 0);

            double s = 1 - xi;
            double r = xi / s;
            double B = wgt * pow(s, m);
            for (uint a2 = 0; a2 <= m; a2++)
            {
                for (uint a1 = 0; a1 <= n; a1++)
                {
                    for (uint k = 0; k < q; k++)
                    {
                        uint index_a1a2k = position_q(a1, a2, k, q);
                        uint index_a1jk = position_q(a1, j, k, q);

                        for (uint ell = 0; ell < nb_Array; ell++)
                            BMomentInter.at(index_a1a2k, ell) += B * Bmoment.at(index_a1jk, ell);
                    }
                }
                B *= r * (m - a2) / (1 + a2);
            }
        }

        Bmoment.zeros();
        // convert third index (l=1)
        for (uint k = 0; k < q; k++)
        {
            double xi = quadraWN.at(k, 1);
            double wgt = quadraWN.at(k, 0);

            double s = 1 - xi;
            double r = xi / s;
            double B = wgt * pow(s, p);
            for (uint a3 = 0; a3 <= p; a3++)
            {
                for (uint a1 = 0; a1 <= n; a1++)
                {
                    for (uint a2 = 0; a2 <= m; a2++)
                    {
                        uint index_a1a2k = position_q(a1, a2, k, q);
                        uint index_a = element.position({a1, a2, a3}, max_nm);

                        for (uint ell = 0; ell < nb_Array; ell++)
                            Bmoment.at(index_a, ell) += B * BMomentInter.at(index_a1a2k, ell);
                    }
                }
                B *= r * (p - a3) / (1 + a3);
            }
        }
    }
    return Bmoment;
}