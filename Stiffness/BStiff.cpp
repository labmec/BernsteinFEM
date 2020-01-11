#include "StiffM.h"

// default constructor
BStiff::BStiff(uint32_t q, uint32_t n)
    : Matrix(), BinomialMat(n + 1, n + 1)
{
    this->q = q;
    this->n = n;
    lenBinomialMat = n + 1;
    computeBinomials();
}

// copy constructor
BStiff::BStiff(const BStiff &cp)
{
    q = cp.q;
    n = cp.n;
    lenBinomialMat = cp.lenBinomialMat;
    BinomialMat = cp.BinomialMat;
    lenStiff = cp.lenStiff;
    Matrix = cp.Matrix;
}

// copy assignment operator (unnecessary)
BStiff &BStiff::operator=(const BStiff &cp)
{
    if (this != &cp)
    {
        q = cp.q;
        n = cp.n;
        lenBinomialMat = cp.lenBinomialMat;
        BinomialMat = cp.BinomialMat;
        lenStiff = cp.lenStiff;
        Matrix = cp.Matrix;
    }
    return *this;
}

void BStiff::computeBinomials()
{
    BinomialMat(0, 0) = BinomialMat(0, 1) = BinomialMat(1, 0) = 1;
    for (uint32_t k = 1; k < lenBinomialMat; k++)
    {
        for (uint32_t l = 1; l < lenBinomialMat; l++)
        {
            BinomialMat(k, l) = BinomialMat(k, l - 1) + BinomialMat(k - 1, l);
        }
    }
}

// returns the length of the matrix
uint32_t BStiff::length()
{
    return lenStiff;
}

uint32_t BStiff::getPOrder()
{
	return n;
}

// returns the value of Matrix[i][j]
double BStiff::getMatrixValue(uint32_t i, uint32_t j)
{
    return Matrix(i, j);
}

// returns the matrix
const TPZFMatrix<REAL> &BStiff::getMatrix()
{
    return Matrix;
}

// other methods

// assign 0 to all elements of the matrix
void BStiff::zero()
{
	Matrix.Zero();
}

// copyBStiff implementations
template<>
BStiff* copyBStiff<Element_t::LinearEl>(BStiff const& copy) {
	return new BStiff1D(dynamic_cast<BStiff1D const&>(copy));
}

template<>
BStiff* copyBStiff<Element_t::TriangularEl>(BStiff const& copy) {
	return new BStiff2DTri(dynamic_cast<BStiff2DTri const&>(copy));
}

template<>
BStiff* copyBStiff<Element_t::QuadrilateralEl>(BStiff const& copy) {
	//return new BStiff2DQuad(dynamic_cast<BStiff2DQuad const&>(copy));
	return nullptr;
}

template<>
BStiff* copyBStiff<Element_t::TetrahedronEl>(BStiff const& copy) {
	//return new BStiff3DTetra(dynamic_cast<BStiff3DTetra const&>(copy));
	return nullptr;
}

template<>
BStiff* copyBStiff<Element_t::CubeEl>(BStiff const& copy) {
	//return new BStiff3DCube(dynamic_cast<BStiff3DTetra const&>(copy));
	return nullptr;
}