#pragma once

#include <vector>
#include <cmath>

namespace VectorUtils
{
    void multiplyVecs(const std::vector<double>& f, const std::vector<double>& g, std::vector<double>& fg, size_t size);

    double weight(int i, int n, double step);

    double jacobian(double x, int dim);

    double getDot(const std::vector<double>& f, const std::vector<double>& g, size_t size, double step, double origin,
                  int dim);

    double getIntNorm(const std::vector<double>& f, size_t size, double step, double origin, int dim);
}