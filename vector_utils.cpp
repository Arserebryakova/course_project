#include "vector_utils.hpp"

void VectorUtils::multiplyVecs(const std::vector<double>& f, const std::vector<double>& g, std::vector<double>& fg,
                               size_t size) {
    for (size_t i = 0; i < size; ++i)
        fg[i] = f[i] * g[i];
}

double VectorUtils::weight(int i, int n, double step) {
    return i == 0 || i == n - 1 ? step / 2 : step;
}

double VectorUtils::jacobian(double x, int dim) {
    return dim == 1 ? 1.0 :
           dim == 2 ? x :
           x * x;
}

double
VectorUtils::getDot(const std::vector<double>& f, const std::vector<double>& g, size_t size, double step, double origin,
                    int dim) {
    double res = 0.0;
    double x = origin;

    for (size_t i = 0; i < size; ++i) {
        res += f[i] * g[i] * weight(i, size, step) * jacobian(x, dim);
        x += step;
    }

    double dimCoeff = dim == 1 ? 2.0 : dim == 2 ? 2.0 * M_PI : 4.0 * M_PI;
    return res * dimCoeff;
}

double VectorUtils::getIntNorm(const std::vector<double>& f, size_t size, double step, double origin, int dim) {
    double res = 0.0;
    double x = origin;

    for (size_t i = 0; i < size; ++i) {
        res += std::abs(f[i]) * weight(i, size, step) * jacobian(x, dim);
        x += step;
    }

    double dimCoeff = dim == 1 ? 2.0 : dim == 2 ? 2.0 * M_PI : 4.0 * M_PI;
    return res * dimCoeff;
}
