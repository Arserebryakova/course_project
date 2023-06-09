#pragma once

#include <math.h>
#include <vector>
#include "problem.hpp"
#include "vector_utils.hpp"

class EulerSolver
{
public:
    Result solve(const Problem& p);

private:
    double N;                     /* first moment */
    std::vector<double> C;        /* second moment samples */

    std::vector<double> m;        /* birth kernel samples */
    std::vector<double> w;        /* death kernel samples */

    std::vector<double> mC;       /* samples of [m * C] */
    std::vector<double> wC;       /* samples of [w * C] */
    std::vector<double> CwC;      /* samples of [Cw * C] */
    std::vector<double> w_mult_C; /* samples of wC */

    /* get zero padded samples of birth and death kernels and for C */
    void getVectors(const Problem& p);

    void getConvolutions(const Problem& p);

    std::vector<double> convolve(const std::vector<double>& f, const std::vector<double>& g, const Problem& p);

    /* convolving function */
};
