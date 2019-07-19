#include "Moments.h"
#include "JacobiGaussNodes.h"

uint position_q(uint i, uint j, uint k, uint q) { return i * q * q + j * q + k; }

REAL det3(TPZFMatrix<REAL> & mat)
{
	// compute determinant by Laplace formula
	REAL det = 0;
	det += mat(0, 0) * (mat(1, 2) * mat(2, 1) - mat(1, 1) * mat(2, 2));
	det += mat(0, 1) * (mat(1, 0) * mat(2, 2) - mat(1, 2) * mat(2, 0));
	det += mat(0, 2) * (mat(1, 1) * mat(2, 0) - mat(1, 0) * mat(2, 1));
	return det;
}

BMoment3DCube::BMoment3DCube(uint q, uint n, const Element<Element_t::CubeEl> &el)
    : BMoment3DCube(q, n, n, n, el) {}

BMoment3DCube::BMoment3DCube(uint q, uint n, uint m, uint p, const Element<Element_t::CubeEl> &el)
    : BMoment(q, n, el), 
    BMomentInter()
{
    uint max = MAX(MAX(n, m), p);
    max++;
    lenMoments = max * max * max; // (max + 1) ^ 3
    lenCval = q * q * q;
    Bmoment.resize(lenMoments);
    BMomentInter.resize(lenMoments);
    Cval.resize(lenCval);
    assignQuadra();
}

// copy constructor
BMoment3DCube::BMoment3DCube(const BMoment3DCube &cp)
    : BMoment(cp.q, cp.n, cp.element), BMomentInter(cp.BMomentInter)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
}

BMoment3DCube &BMoment3DCube::operator=(const BMoment3DCube &cp)
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
void BMoment3DCube::assignQuadra()
{
    for (uint k = 0; k < q; k++)
    {
        intPoints[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        intWeights[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment3DCube::loadFunctionDef()
{
    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; i++)
        {
            for (uint k = 0; k < q; k++)
            {
                TPZFMatrix<REAL> jac(3, 3);
                TPZFMatrix<REAL> X = element.mapToElement({intPoints[i], intPoints[j], intPoints[k]}, jac);
                REAL dX = det3(jac);
                uint index_ijk = position_q(i, j, k, q);
                Cval[index_ijk] = f(X[0], X[1], X[2]) * dX;
            }
        }
    }
}

inline
TPZFMatrix<REAL> BMoment3DCube::getIntegrationPoints()
{
    TPZFMatrix<REAL> points(lenCval, 3);

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            for (uint k = 0; k < q; k++)
            {
                TPZFMatrix<REAL> X = element.mapToElement({intPoints[i], intPoints[j], intPoints[k]});
                uint index_ijk = position_q(i, j, k, q);
                points(index_ijk, 0) = X[0];
                points(index_ijk, 1) = X[1];
                points(index_ijk, 2) = X[2];
            }
        }
    }

    return points;
}

TPZVec<REAL> &BMoment3DCube::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        if (!functVal)
            loadFunctionDef();

        uint max_nm = MAX(MAX(n, m), p);
        uint max_nq = MAX(max_nm, q);
		Bmoment.Fill(0.0);
        BMomentInter.Fill(0.0);

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
                    for (uint k = 0; j < q; k++)
                    {
                        uint index_a1jk = position_q(a1, j, k, max_nq);
                        uint index_ijk = position_q(i, j, k, q);

						Bmoment[index_a1jk] += B * Cval[index_ijk];
                    }
                }
                B *= r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=2)
        for (uint j = 0; j < q; j++)
        {
			double xi = intPoints[j];
			double wgt = intWeights[j];

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

						BMomentInter[index_a1a2k] += B * Bmoment[index_a1jk];
                    }
                }
                B *= r * (m - a2) / (1 + a2);
            }
        }

		Bmoment.Fill(0.0);
        // convert third index (l=1)
        for (uint k = 0; k < q; k++)
        {
			double xi = intPoints[k];
			double wgt = intWeights[k];

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
                        uint index_a = element.position({a1, a2, a3});

						Bmoment[index_a] += B * BMomentInter[index_a1a2k];
                    }
                }
                B *= r * (p - a3) / (1 + a3);
            }
        }
    }
    return Bmoment;
}