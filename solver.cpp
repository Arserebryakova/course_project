#include "solver.hpp"
#include <functional>
#include <numeric>
#include <iostream>

inline std::vector<double> calculateWeightedJacobians(const Problem& p) {
    double dimCoeff = p.dimension() == 1 ? 2.0 : p.dimension() == 2 ? 2.0 * M_PI : 4.0 * M_PI;
    std::vector<double> res(p.nodes(), p.step() * dimCoeff);
    res[0] = res[p.nodes() - 1] = p.step() / 2;
    double x = p.origin();
    for (size_t i = 0; i < p.nodes(); ++i) {
        if (p.dimension() == 2) res[i] *= x;
        else res[i] *= x * x;
        x += p.step();
    }
    return res;
}

Result EulerSolver::solve(const Problem& p) {
    auto vecs = Samples(p.nodes());
    auto weightedJacobians = calculateWeightedJacobians(p);
    precalcSamples(vecs, p, weightedJacobians);

    for (size_t i = 0; i < p.iters(); ++i) {
        std::transform(vecs.w.begin(), vecs.w.end(), vecs.C.begin(), vecs.w_mult_C.begin(), std::multiplies<double>());
        double weightedDotProduct = std::inner_product(vecs.w_mult_C.begin(), vecs.w_mult_C.end(),
                                                       weightedJacobians.begin(), 0);
        vecs.N = (p.b() - p.d()) / (weightedDotProduct + p.s());
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

void EulerSolver::precalcSamples(Samples& vecs, const Problem& p, const std::vector<double>& weightedJacobians) {
    double x = p.origin();

    for (size_t i = 0; i < p.nodes(); ++i) {
        vecs.m[i] = p.getKernels().m(x);
        vecs.w[i] = p.getKernels().w(x);
        x += p.step();
    }

    auto absMultiplication = [](const auto& a, const auto& b) { return std::abs(a * b); };
    double nm = std::inner_product(vecs.m.begin(), vecs.m.end(), weightedJacobians.begin(), 0, std::plus<double>(), absMultiplication); // weighted norm of m
    double nw = std::inner_product(vecs.w.begin(), vecs.w.end(), weightedJacobians.begin(), 0, std::plus<double>(), absMultiplication); // weighted norm of w

    for (size_t i = 0; i < p.nodes(); ++i) {
        vecs.m[i] *= p.b() / nm;
        vecs.C[i] = vecs.w[i] = vecs.w[i] * p.s() / nw;
    }
}

void EulerSolver::getConvolutions(Samples& vecs, const Problem& p) {
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