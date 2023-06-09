#include "solver.hpp"
#include <iostream>
#include <iomanip>
#include <fstream>

Result EulerSolver::solve(const Problem& p)
{
    Result res;
    getVectors(p);

    for (size_t i = 0; i < p.iters(); ++i)
    {
        N = (p.b() - p.d()) / (VectorUtils::getDot(C, w, p.nodes(), p.step(), p.origin(), p.dimension()) + p.s());
        getConvolutions(p);

        for (size_t j = 0; j < p.nodes(); ++j)
        {
            C[j] = (m[j] / N - w[j] + mC[j] + p.b() - p.d() -
                    N / (p.alpha() + p.gamma()) *
                    (p.alpha() * (p.b() - p.d()) / N +
                     p.beta() * (wC[j] + CwC[j]) +
                     p.gamma() * ((p.b() - p.d()) / N + wC[j] + CwC[j]))) /
                   (p.d() + w[j] + N / (p.alpha() + p.gamma()) *
                                   p.alpha() * (p.b() - p.d()) / N + p.beta() * p.s());
        }
        std::cout << i << " " << N << " " << N * N * (C[0] + 1) << " " << std::endl;
    }

    /* correcting second moment */
    double NN = N * N;
    for (size_t i = 0; i < p.nodes(); ++i)
    {
        C[i] = NN * (C[i] + 1);
    }

    res.N = N;
    res.C.assign(p.nodes(), 0);
    res.dim = p.dimension();
    res.n_count = p.nodes();
    VectorUtils::copy(res.C, C, p.nodes());

    return res;
}

void EulerSolver::getVectors(const Problem& p)
{
    w.assign(p.nodes(), 0);
    m.assign(p.nodes(), 0);
    C.assign(p.nodes(), 0);

    mC.assign(p.nodes(), 0);
    wC.assign(p.nodes(), 0);
    CwC.assign(p.nodes(), 0);
    w_mult_C.assign(p.nodes(), 0);

    double x = p.origin();

    for (size_t i = 0; i < p.nodes(); ++i)
    {
        m[i] = p.getKernels().m(x);
        w[i] = p.getKernels().w(x);
        x += p.step();
    }

    double nm = VectorUtils::getIntNorm(m, p.nodes(), p.step(), p.origin(), p.dimension());
    double nw = VectorUtils::getIntNorm(w, p.nodes(), p.step(), p.origin(), p.dimension());

    for (size_t i = 0; i < p.nodes(); ++i)
    {
        m[i] *= p.b() / nm;
        C[i] = w[i] = w[i] * p.s() / nw;
    }
}

void EulerSolver::getConvolutions(const Problem& p)
{
    VectorUtils::multiplyVecs(C, w, w_mult_C, p.nodes());

    if (p.dimension() == 3)
    {
        double x = p.origin();

        for (size_t i = 0; i < p.nodes(); ++i)
        {
            w_mult_C[i] *= 4 * M_PI * std::abs(x);
            x += p.step();
        }
    }
    mC = convolve(m, C, p);
    wC = convolve(w, C, p);
    CwC = convolve(C, w_mult_C, p);
}

std::vector<double> EulerSolver::convolve(const std::vector<double>& f, const std::vector<double>& g, const Problem& p)
{
    std::vector<double> fg(f.size(), 0);
    for (size_t i = 0; i < p.nodes(); ++i)
    {
        for (size_t j = 0; j < i; ++j)
            fg[i] += f[i - j] * g[j] * p.step();
    }
    return fg;
}
