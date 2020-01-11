#include "StiffM.h"

BStiff1D::BStiff1D(uint32_t q, uint32_t n, const Element<Element_t::LinearEl> &element)
    : BStiff(q, n), BMoment1D(q, 2 * (n - 1), element)
{
    lenStiff = n + 1;
	Matrix.Resize(lenStiff, lenStiff);
    this->element.setPermutationPOrder(n);
}

BStiff1D::BStiff1D(const BStiff1D& cp) :
	BStiff(cp),	BMoment1D(cp) {
	lenStiff = cp.lenStiff;
	Matrix = cp.Matrix;
	element = cp.element;
}

BStiff1D& BStiff1D::operator=(const BStiff1D& cp) {
	if (this != &cp) {
		BStiff::operator=(cp);
		BMoment1D::operator=(cp);
	}
	return *this;
}

void BStiff1D::computeMatrix()
{
    computeMoments();
    computeBinomials();

    const auto n = BStiff::n;
    const auto a = element.getVertices()(0, 0), b = element.getVertices()(0, 1);

    const auto dist = b - a;
    const auto Const = 1. / BinomialMat(n - 1, n - 1) / pow(dist, n);
    const auto n2 = n * n;

    for (uint i = 0; i < lenStiff - 1; i++)
    {
        for (uint j = 0; j < lenStiff - 1; j++)
        {
            double w = BinomialMat(i, j) * Const;
            w *= BinomialMat(n - i - 1, n - j - 1);

        	w *= n2 * dist * Bmoment[i + j];
            const auto _i = element.position({i});
            const auto _j = element.position({j});
            const auto I = element.position({i + 1});
            const auto J = element.position({j + 1});

            Matrix(_i, _j) += w;
            Matrix(I, J) += w;
            Matrix(I, _j) -= w;
            Matrix(_i, J) -= w;
        }
    }
}