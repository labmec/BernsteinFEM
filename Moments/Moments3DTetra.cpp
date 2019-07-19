#include "Moments.h"
#include "JacobiGaussNodes.h"

// helps indexing quadrature points vectors
int position_q(uint i, uint j, uint k, uint q) { return i * q * q + j * q + k; }

REAL Volume(TPZFMatrix<REAL> const &v);

BMoment3DTetra::BMoment3DTetra(uint q, uint n, const Element<Element_t::TetrahedronEl> &element)
    : BMoment(q, n, element), BMomentInter()
{
    int max = MAX(n + 1, q);
    n++;
    lenMoments = n * n * n;  // FIXME: it's possible to do this with (n + 1) * (n + 2) * (n + 3) / 6
    lenCval = q * q * q;
    Bmoment.resize(lenMoments); // FIXME: initial size for this matrix might need to bigger
    BMomentInter.resize(max * max * max);
    Cval.resize(lenCval);
    assignQuadra();
}

BMoment3DTetra::BMoment3DTetra(const BMoment3DTetra &cp)
    : BMoment(cp.q, cp.n, cp.element),
      BMomentInter(cp.BMomentInter)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
}

BMoment3DTetra &BMoment3DTetra::operator=(const BMoment3DTetra &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        lenMoments = cp.lenMoments;
        lenCval = cp.lenCval;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
    }
    return *this;
}

inline
void BMoment3DTetra::assignQuadra()
{
    for (uint k = 0; k < q; k++)
    {
        intPoints[k] = (1.0 + jacobi2_xi(q, k)) * 0.5;
        intWeights[k] = jacobi2_w(q, k) * 0.25;
    }

    for (uint k = q; k < 2 * q; k++)
    {
        intPoints[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        intWeights[k] = legendre_w(q, k) * 0.5;
    }

    for (uint k = 2 * q; k < 3 * q; k++)
    {
        intPoints[k] = (1.0 + jacobi_xi(q, k)) * 0.5;
        intWeights[k] = jacobi_w(q, k) * 0.25;
    }
}

inline
void BMoment3DTetra::loadFunctionDef()
{
    TPZFMatrix<REAL> points(getIntegrationPoints());
    Cval.resize(lenCval);
    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            for (uint k = 0; k < q; k++)
            {
                uint index_ijk = position_q(i, j, k, q);
                Cval[index_ijk] = f(points(i, 0), points(j, 1), points(k, 2));
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
TPZFMatrix<REAL> BMoment3DTetra::getIntegrationPoints()
{
    TPZFMatrix<REAL> points(q * q * q , 3); // vector with q * q  elements

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            for (uint k = 0; k < q; k++)
            {
                TPZFMatrix<REAL> v = element.mapToElement({intPoints[i], intPoints[j + q], intPoints[k + q * 2]});
                uint index_ijk = position_q(i, j, k, q);
                points(index_ijk, 0) = v[0];
                points(index_ijk, 1) = v[1];
                points(index_ijk, 2) = v[2];
            }
        }
    }

    return points;
}

// compute the b-moments
TPZVec<REAL> &BMoment3DTetra::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'computeMoments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'computeMoments()\'\n";
    else
    {
        if (!functVal)
            loadFunctionDef();

		Bmoment.Fill(0.0);
        BMomentInter.Fill(0.0);

        uint m = MAX(n, q - 1);

        double scalingConst = 6 * Volume(element.getVertices()); // FIXME: check this value out

        // convert first index (l=3)
        for (uint i = 0; i < q; i++)
        {
            double xi = intPoints[i];
			double wgt = intWeights[i];

            double s = 1 - xi;
            double r = xi / s;
            double B = wgt * pow(s, n);
            for (uint a1 = 0; a1 <= n; a1++)
            {
                for (uint j = 0; j < q; j++)
                {
                    for (uint k = 0; k < q; k++)
                    {
                        uint index_a1jk = position_q(a1, j, k, m);
                        uint index_ijk = position_q(i, j, k, q);

                         Bmoment[index_a1jk] += scalingConst * B * Cval[index_ijk];
					     // Bmoment(a1 * (n + 1) * (n + 1) + j * q + k) 
                         // += scalingConst * B *Cval(i * q * q + j * q + k, ell);
                    }
                }
                B *= r * (n - a1) / (1 + a1);
            }
        }

        // convert second index
        for (uint j = 0; j < q; j++)
        {
			double xi = intPoints[j + q];
			double wgt = intWeights[j + q];

            double s = 1 - xi;
            double r = xi / s;
            for (uint a1 = 0; a1 <= n; a1++)
            {
                double B = wgt * pow(s, n - a1);
                for (uint a2 = 0; a2 <= n - a1; a2++)
                {
                    for (uint k = 0; k < q; k++)
                    {
                        uint index_a1jk = position_q(a1, j, k, m);
                        uint index_a1a2k = position_q(a1, a2, k, m);

						for (uint ell = 0; ell; ell++)
							BMomentInter[index_a1a2k] += B * Bmoment[index_a1jk];
                    }
                    B *= r * (n - a1 - a2) / (1 + a2);
                }
            }
        }

		Bmoment.Fill(0.0);
        // conver third index
        for (uint k = 0; k < q; k++)
        {
			double xi = intPoints[k + q * 2];
			double wgt = intWeights[k + q * 2];

            double s = 1 - xi;
            double r = xi / s;
            for (uint a1 = 0; a1 <= n; a1++)
            {
                for (uint a2 = 0; a2 <= n - a1; a2++)
                {
                    double B = wgt * pow(s, n - a1 - a2);
                    for (uint a3 = 0; a3 <= n - a1 - a2; a3++)
                    {
                        uint index_a1a2k = position_q(a1, a2, k, m);
                        uint index = element.position({a1, a2, a3});

                        for (uint ell = 0; ell ; ell++)
                            Bmoment[index] += B * BMomentInter[index_a1a2k];

                        B *= r * (n - a1 - a2 - a3) / (1 + a3);
                    }
                }
            }
        }
    }
    return Bmoment;
}

// computes volume of tetrahedron defined by vertices v
REAL Volume(TPZFMatrix<REAL> const &v)
{
	return 0;
    // TODO: calc new volume
}