#include "solver.hpp"
#include <iostream>

Result EulerSolver::solve(const Problem& p) {
    auto vecs = Samples(p.nodes());
    precalcSamples(vecs, p);

    for (size_t i = 0; i < p.iters(); ++i) {
        vecs.N = (p.b() - p.d()) /
                 (VectorUtils::getDot(vecs.C, vecs.w, p.nodes(), p.step(), p.origin(), p.dimension()) + p.s());
        getConvolutions(vecs, p);

        for (size_t j = 0; j < p.nodes(); ++j) {
            vecs.C[j] = (vecs.m[j] / vecs.N - vecs.w[j] + vecs.mC[j] + p.b() - p.d() -
                         vecs.N / (p.alpha() + p.gamma()) *
                         (p.alpha() * (p.b() - p.d()) / vecs.N +
                          p.beta() * (vecs.wC[j] + vecs.CwC[j]) +
                          p.gamma() * ((p.b() - p.d()) / vecs.N + vecs.wC[j] + vecs.CwC[j]))) /
                        (p.d() + vecs.w[j] + vecs.N / (p.alpha() + p.gamma()) *
                                             p.alpha() * (p.b() - p.d()) / vecs.N + p.beta() * p.s());
        }
        std::cout << i << " " << vecs.N << " " << vecs.N * vecs.N * (vecs.C[0] + 1) << " " << std::endl;
    }

    /* correcting second moment */
    double NN = vecs.N * vecs.N;
    for (size_t i = 0; i < p.nodes(); ++i)
        vecs.C[i] = NN * (vecs.C[i] + 1);

    return {vecs.N, p.dimension(), p.nodes(), vecs.C};
}

void EulerSolver::precalcSamples(Samples& vecs, const Problem& p) {
    double x = p.origin();

    for (size_t i = 0; i < p.nodes(); ++i) {
        vecs.m[i] = p.getKernels().m(x);
        vecs.w[i] = p.getKernels().w(x);
        x += p.step();
    }

    double nm = VectorUtils::getIntNorm(vecs.m, p.nodes(), p.step(), p.origin(), p.dimension());
    double nw = VectorUtils::getIntNorm(vecs.w, p.nodes(), p.step(), p.origin(), p.dimension());

    for (size_t i = 0; i < p.nodes(); ++i) {
        vecs.m[i] *= p.b() / nm;
        vecs.C[i] = vecs.w[i] = vecs.w[i] * p.s() / nw;
    }
}

void EulerSolver::getConvolutions(Samples& vecs, const Problem& p) {
    VectorUtils::multiplyVecs(vecs.C, vecs.w, vecs.w_mult_C, p.nodes());

    if (p.dimension() == 3) {
        double x = p.origin();

        for (size_t i = 0; i < p.nodes(); ++i) {
            vecs.w_mult_C[i] *= 4 * M_PI * std::abs(x);
            x += p.step();
        }
    }
    vecs.mC = convolve(vecs.m, vecs.C, p);
    vecs.wC = convolve(vecs.w, vecs.C, p);
    vecs.CwC = convolve(vecs.C, vecs.w_mult_C, p);
}

std::vector<double>
EulerSolver::convolve(const std::vector<double>& f, const std::vector<double>& g, const Problem& p) {
    std::vector<double> fg(f.size(), 0);
    for (size_t i = 0; i < p.nodes(); ++i) {
        for (size_t j = 0; j < i; ++j)
            fg[i] += f[i - j] * g[j] * p.step();
    }
    return fg;
}

Samples::Samples(size_t n_nodes) : C(n_nodes), m(n_nodes), w(n_nodes), mC(n_nodes), wC(n_nodes), CwC(n_nodes),
                                   w_mult_C(n_nodes) {

}