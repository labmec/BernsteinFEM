#include "Evaluation.h"
#include "JacobiGaussNodes.h"

BEval2DTri::BEval2DTri(uint q, uint n, Element<Element_t::TriangularEl> const &el)
    : BEval(q, n, el), eval_inter(q * q)
{
    bbvec_len = (n + 1) * (n + 1);
    eval.set_size(q * q);
}

BEval2DTri::BEval2DTri(uint q, uint n, arma::vec const &cVec, Element<Element_t::TriangularEl> const &el)
    : BEval(q, n, cVec, el), eval_inter(q * q)
{
    bbvec_len = (n + 1) * (n + 1);
    eval.set_size(q * q);
}

arma::vec &BEval2DTri::computeEvaluation()
{
    eval.zeros();
    eval_inter.zeros();

    // convert first index
    for (uint i2 = 0; i2 < q; i2++)
    {
        double xi = (1.0 + legendre_xi(q, i2)) * 0.5;
        double s = 1 - xi;
        double r = xi / s;

        for (uint a1 = 0; a1 <= n; a1++)
        {
            double w = pow(s, n - a1);
            for (uint a2 = 0; a2 <= n - a1; a2++)
            {
                eval_inter.at(i2 * q + a1) += w * BBVec.at(element.position({a1, a2}));
                w *= r * (n - a1 - a2) / (1 + a2);
            }
        }
    }
    
    // convert second index
    for (uint i1 = 0; i1 < q; i1++)
    {
        double xi = (1.0 + jacobi_xi(q, i1)) * 0.5;
        double s = 1 - xi;
        double r = xi / s;

        double w = pow(s, n);
        for (uint a1 = 0; a1 <= n; a1++)
        {
            for (uint i2 = 0; i2 < q; i2++)
            {
                eval.at(i1 * q + i2) += eval_inter.at(i2 * q + a1);
            }
            w *= r * (n - a1) / (1 + a1);
        }
    }
    evaluated = true;
    return eval;
}

double BEval2DTri::L2Norm()
{
    if (!evaluated)
        computeEvaluation();

    double norm = 0;

    for (uint i = 0; i < q; i++)
    {
        double w1 = legendre_w(q, i) * 0.5;
        for (uint j = 0; j < q; j++)
        {
            double w2 = w1 * (jacobi_w(q, j) * 0.5);
            double f = eval.at(i * q + j);
            norm += w2 * f * f;
        }
    }

    return sqrt(norm);
}