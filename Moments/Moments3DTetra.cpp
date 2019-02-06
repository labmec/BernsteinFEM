#include "Moments.h"
#include "JacobiGaussNodes.h"

// helps indexing quadrature points vectors
int position_q(int i, int j, int k, int q) { return i * q * q + j * q + k; }

double Volume(arma::mat const &v);

BMoment3DTetra::BMoment3DTetra(int q, int n, const Element<Element_t::TetrahedronEl> &element, int nb_Array)
    : BMoment(q, n, element, nb_Array), BMomentInter()
{
    int max = MAX(n + 1, q);
    n++;
    lenMoments = n * n * n;  // FIXME: it's possible to do this with (n + 1) * (n + 2) * (n + 3) / 6
    lenCval = q * q * q;
    Bmoment.set_size(lenMoments, nb_Array); // FIXME: initial size for this matrix might need to bigger
    BMomentInter.set_size(max * max * max, nb_Array);
    Cval.set_size(lenCval, nb_Array);
    assignQuadra();
}

BMoment3DTetra::BMoment3DTetra(const BMoment3DTetra &cp)
    : BMoment(cp.q, cp.n, cp.element, cp.nb_Array),
      BMomentInter(cp.BMomentInter)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
    quadraWN = cp.quadraWN;
}

BMoment3DTetra &BMoment3DTetra::operator=(const BMoment3DTetra &cp)
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
void BMoment3DTetra::assignQuadra()
{
    quadraWN.zeros(q, 6);
    double *x = quadraWN.colptr(1);
    double *w = quadraWN.colptr(0);

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + jacobi2_xi(q, k)) * 0.5;
        w[k] = jacobi2_w(q, k) * 0.25;
    }

    x = quadraWN.colptr(3);
    w = quadraWN.colptr(2);

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        w[k] = legendre_w(q, k) * 0.5;
    }

    x = quadraWN.colptr(5);
    w = quadraWN.colptr(4);

    for (int k = 0; k < q; k++)
    {
        x[k] = (1.0 + jacobi_xi(q, k)) * 0.5;
        w[k] = jacobi_w(q, k) * 0.25;
    }
}

inline
void BMoment3DTetra::loadFunctionDef()
{
    arma::mat points(getIntegrationPoints());
    Cval.set_size(lenCval, nb_Array);
    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < q; k++)
            {
                int index_ijk = position_q(i, j, k, q);
                Cval.at(index_ijk, 0) = f(points.at(i, 0), points.at(j, 1), points.at(k, 2));
            }
        }
    }
    fValSet = true;
}

// returns the vectors with the integration points (x, y) over the object's element, following the moments organization
// usage: points = getIntegrationPoints(); then
// points(i, 0) == x i-th coordinate
// points(i, 1) == y i-th coordinate
// points(i, 2) == z i-th coordinate
inline
arma::mat BMoment3DTetra::getIntegrationPoints()
{
    arma::mat points(q * q * q * nb_Array, 3); // vector with q * q * nb_Array elements

    for (int i = 0; i < q; i++)
    {
        for (int j = 0; j < q; j++)
        {
            for (int k = 0; k < q; k++)
            {
                arma::mat v = element.mapToElement({quadraWN.at(i, 1), quadraWN.at(j, 3), quadraWN(k, 5)});
                int index_ijk = position_q(i, j, k, q);
                points(index_ijk, 0) = v[0];
                points(index_ijk, 1) = v[1];
                points(index_ijk, 2) = v[2];
            }
        }
    }

    return points;
}

// compute the b-moments
arma::mat &BMoment3DTetra::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'computeMoments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'computeMoments()\'\n";
    else
    {
        if (!functVal)
            loadFunctionDef();

        Bmoment.zeros();
        BMomentInter.zeros();

        int m = MAX(n, q - 1);

        double scalingConst = 6 * Volume(element.getVertices()); // FIXME: check this value out

        // convert first index (l=3)
        for (int i = 0; i < q; i++)
        {
            double xi = quadraWN.at(i, 1);
            double wgt = quadraWN.at(i, 0);

            double s = 1 - xi;
            double r = xi / s;
            double B = wgt * pow(s, n);
            for (int a1 = 0; a1 <= n; a1++)
            {
                for (int j = 0; j < q; j++)
                {
                    for (int k = 0; k < q; k++)
                    {
                        int index_a1jk = position_q(a1, j, k, m);
                        int index_ijk = position_q(i, j, k, q);

                        for (int ell = 0; ell < nb_Array; ell++)
                            Bmoment(index_a1jk) += scalingConst * B * Cval(index_ijk);
                            // Bmoment(a1 * (n + 1) * (n + 1) + j * q + k) 
                            // += scalingConst * B *Cval(i * q * q + j * q + k, ell);
                    }
                }
                B *= r * (n - a1) / (1 + a1);
            }
        }

        // convert second index
        for (int j = 0; j < q; j++)
        {
            double xi = quadraWN.at(j, 3);
            double wgt = quadraWN.at(j, 2);

            double s = 1 - xi;
            double r = xi / s;
            for (int a1 = 0; a1 <= n; a1++)
            {
                double B = wgt * pow(s, n - a1);
                for (int a2 = 0; a2 <= n - a1; a2++)
                {
                    for (int k = 0; k < q; k++)
                    {
                        int index_a1jk = position_q(a1, j, k, m);
                        int index_a1a2k = position_q(a1, a2, k, m);

                        for (int ell = 0; ell < nb_Array; ell++)
                            BMomentInter.at(index_a1a2k, nb_Array) += B * Bmoment.at(index_a1jk, nb_Array);
                    }
                    B *= r * (n - a1 - a2) / (1 + a2);
                }
            }
        }

        Bmoment.zeros();
        // conver third index
        for (int k = 0; k < q; k++)
        {
            double xi = quadraWN.at(k, 5);
            double wgt = quadraWN.at(k, 4);

            double s = 1 - xi;
            double r = xi / s;
            for (int a1 = 0; a1 <= n; a1++)
            {
                for (int a2 = 0; a2 <= n - a1; a2++)
                {
                    double B = wgt * pow(s, n - a1 - a2);
                    for (int a3 = 0; a3 <= n - a1 - a2; a3++)
                    {
                        int index_a1a2k = position_q(a1, a2, k, m);
                        int index = position(a1, a2, a3, n);

                        for (int ell = 0; ell < nb_Array; ell++)
                            Bmoment.at(index, nb_Array) += B * BMomentInter.at(index_a1a2k, nb_Array);

                        B *= r * (n - a1 - a2 - a3) / (1 + a3);
                    }
                }
            }
        }
    }
    return Bmoment;
}

// computes volume of tetrahedron defined by vertices v
double Volume(arma::mat const &v)
{
    arma::mat m(4, 1, arma::fill::ones);
    m.insert_cols(1, v);

    return fabs(arma::det(m)) / 6.;
}