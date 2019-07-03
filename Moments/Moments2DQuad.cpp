#include "Moments.h"
#include "JacobiGaussNodes.h"

using std::vector;

// helps indexing quadrature points vectors
uint position_q(uint i, uint j, int q) { return i * q + j; }

BMoment2DQuad::BMoment2DQuad(uint q, uint n, const Element<Element_t::QuadrilateralEl> &element)
    : BMoment2DQuad(q, n, n, element) {}

BMoment2DQuad::BMoment2DQuad(uint q, uint n, uint m, const Element<Element_t::QuadrilateralEl> &element)
	: BMoment(q, MAX(n, m), element),
	BMomentInter((MAX(MAX(n, m) + 1, q)) * (MAX(MAX(n, m) + 1, q)))
{
    this->m = m < n ? m : n;
    int max = MAX(n, m);
    lenMoments = (max + 1) * (max + 1);
    lenCval = q * q;
    Bmoment.resize(lenMoments);
    Cval.resize(lenCval);
	intPoints.resize(q);
	intWeights.resize(q);
    assignQuadra();
}

// copy constructor
BMoment2DQuad::BMoment2DQuad(const BMoment2DQuad &cp)
    : BMoment(cp.q, cp.n, cp.element),
      BMomentInter(cp.BMomentInter)
{
    lenMoments = cp.lenMoments;
    lenCval = cp.lenCval;
    Bmoment = cp.Bmoment;
    Cval = cp.Cval;
	intPoints = cp.intPoints;
	intWeights = cp.intWeights;
}

// copy assignment operator
BMoment2DQuad &BMoment2DQuad::operator=(const BMoment2DQuad &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        lenMoments = cp.lenMoments;
        lenCval = cp.lenCval;
        Bmoment = cp.Bmoment;
        Cval = cp.Cval;
		intPoints = cp.intPoints;
		intWeights = cp.intWeights;
    }
    return *this;
}

inline
void BMoment2DQuad::assignQuadra()
{
    for (uint k = 0; k < q; k++)
    {
        intPoints[k] = (1.0 + legendre_xi(q, k)) * 0.5;
        intWeights[k] = legendre_w(q, k) * 0.5;
    }
}

inline
void BMoment2DQuad::loadFunctionDef()
{

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            TPZFMatrix<REAL> jac(2, 2);
            TPZFMatrix<REAL> X = element.mapToElement({intPoints[i], intPoints[j]}, jac);
			double dX;
			//jac.DeterminantInverse(dX, jac); // just to get the determinant
			dX = jac(0, 0) * jac(1, 1) - jac(0, 1) * jac(1, 0);
            uint index_ij = position_q(i, j, q);
            Cval[index_ij]= f(X[0], X[1]) * dX;
        }
    }

    fValSet = true;
}

inline
TPZFMatrix<REAL> BMoment2DQuad::getIntegrationPoints()
{
    TPZFMatrix<REAL> points(lenCval, 2);

    for (uint i = 0; i < q; i++)
    {
        for (uint j = 0; j < q; j++)
        {
            TPZFMatrix<REAL> X = element.mapToElement({intPoints[i], intPoints[j]});
            uint index_ij = position_q(i, j, q);
            points(index_ij, 0) = X[0];
            points(index_ij, 1) = X[1];
        }
    }

    return points;
}

TPZVec<REAL> &BMoment2DQuad::computeMoments()
{
    if (functVal && !fValSet)
        std::cerr << "missing function values for computation of the moments in \'compute_moments()\'\n";
    else if (!functVal && !fDefSet)
        std::cerr << "missing function definition for computation of the moments in \'compute_moments()\'\n";
    else
    {
        if (!functVal)
            loadFunctionDef();

        uint max_nm = MAX(n + 1, m + 1);
        uint max_nq = MAX(max_nm, q);
        

		std::fill(Bmoment.begin(), Bmoment.end(), 0.0);
		std::fill(BMomentInter.begin(), BMomentInter.end(), 0.0);

        // convert first index (l=2)
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
                    uint index_a1j = position_q(a1, j, max_nq);

                    uint index_ij = position_q(i, j, q);

					BMomentInter[index_a1j] += B * Cval[index_ij];
                }
                B = B * r * (n - a1) / (1 + a1);
            }
        }

        // convert second index (l=1)
        for (uint i = 0; i < q; i++)
        {
			double xi = intPoints[i];
			double wgt = intWeights[i];

            double s = 1 - xi;
            double r = xi / s;

            double B = wgt * pow(s, m);
            for (uint a2 = 0; a2 <= m; a2++)
            {
                for (uint a1 = 0; a1 <= n; a1++)
                {
                    uint index_a1i = position_q(a1, i, max_nq);
                    uint index_a1a2 = element.position({a1, a2});

					Bmoment[index_a1a2] += B * BMomentInter[index_a1i];
                }
                B = B * r * (m - a2) / (1 + a2);
            }
        }
    }
    return Bmoment;
}
