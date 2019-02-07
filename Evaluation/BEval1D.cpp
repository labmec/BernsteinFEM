#include "Evaluation.h"
#include "JacobiGaussNodes.h"

BEval1D::BEval1D(uint q, uint n, Element<Element_t::LinearEl> const &el)
    : BEval(q, n, el) 
{
    bbvec_len = (n + 1);
    eval.set_size(q);
}

BEval1D::BEval1D(uint q, uint n, arma::vec const &cVec, Element<Element_t::LinearEl> const &el)
    : BEval(q, n, cVec, el)
{
    bbvec_len = (n + 1);
    eval.set_size(q);
}

arma::vec &BEval1D::computeEvaluation()
{
    // set the evaluation vector to 0's
    eval.zeros();

    for (uint i = 0; i < q; i++)
    {
        double xi = (1.0 + legendre_xi(q, i)) * 0.5;
        double s = 1 - xi;
        double r = xi / s;
        double w = pow(s, n);
        for (uint a1 = 0; a1 <= n; a1++)
        {
            eval.at(i) += w * BBVec.at(a1);
            w *= r * (n - a1) / (1 + a1);
        }
    }

    return eval;
}