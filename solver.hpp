#pragma once

#include <math.h>
#include <vector>
#include "problem.hpp"

/* get zero padded samples of birth and death kernels and for C */
struct Samples {
    double N;                     /* first moment */
    std::vector<double> C;        /* second moment samples */

    std::vector<double> m;        /* birth kernel samples */
    std::vector<double> w;        /* death kernel samples */

    std::vector<double> mC;       /* samples of [m * C] */
    std::vector<double> wC;       /* samples of [w * C] */
    std::vector<double> CwC;      /* samples of [Cw * C] */
    std::vector<double> w_mult_C; /* samples of wC */
    Samples(size_t n_nodes);
};

class EulerSolver {
public:
    Result solve(const Problem& p);
private:
    void precalcSamples(Samples& vecs, const Problem& p, const std::vector<double>& weightedJacobians);
    void getConvolutions(Samples& vecs, const Problem& p);
    static std::vector<double> convolve(const std::vector<double>& f, const std::vector<double>& g, const Problem& p);
    /* convolving function */
};
