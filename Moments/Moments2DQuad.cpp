#include "Moments.h"
#include "JacobiGaussNodes.h"

BMoment2DQuad::BMoment2DQuad(int q, int n, Element<Element_t::QuadrilateralEl> element, int nb_Array)
    : BMoment(q, n, element, nb_Array),
      BMomentInter((MAX(n + 1, q)) * (MAX(n + 1, q)), nb_Array, arma::fill::zeros)
{
    lenMoments = (n + 1) * (n + 1);
    lenCval = q * q;
    Bmoment.set_size(lenMoments, nb_Array);
    quadraWN.set_size(lenCval);
    Cval.set_size(lenCval, nb_Array);
    assignQuadra();
}

inline
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

inline
void BMoment2DQuad::loadFunctionDef()
{

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            arma::mat jac(2, 2, arma::fill::none);
            arma::mat X = element.mapToElement({quadraWN.at(i, 1), quadraWN.at(j, 1)}, jac);
            double dX = det(jac);
            int index_ij = position_q(i, j, q);
            Cval.at(index_ij, 0) = f(X[0], X[1]) * dX;
        }
    }

    fValSet = true;
}

inline
arma::mat BMoment2DQuad::getIntegrationPoints()
{
    arma::mat points(lenCval, 2, arma::fill::none);

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            arma::mat X = element.mapToElement({quadraWN.at(i, 1), quadraWN.at(j, 1)});
            int index_ij = position_q(i, j, q);
            points.at(index_ij, 0) = X[0];
            points.at(index_ij, 1) = X[1];
        }
    }

    return points;
}

inline
void BMoment2DQuad::computeMoments()
{
    if (functVal == 1 && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (functVal == 0 && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        if (functVal == 0)
            loadFunctionDef();

        int max_nq = MAX(n, q - 1);
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
                    int index_a1j = position_q(a1, j, max_nq);

                    int index_ij = position_q(i, j, q);

                    for (int ell = 0; ell < nb_Array; ell++)
                        BMomentInter(index_a1j, ell) += B * Cval(index_ij, ell);
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
                double B = wgt * pow(s, n);
                int index_a1i = position_q(a1, i, max_nq);
                for (int a2 = 0; a2 <= n; a2++)
                {
                    int index_a1a2 = position(a1, a2, n);

                    for (int ell = 0; ell < nb_Array; ell++)
                        Bmoment(index_a1a2, 0) += B * BMomentInter(index_a1i, 0);
                    
                    B = B * r * (n - a2) / (1 + a2);
                }
            }
        }

    }
}
