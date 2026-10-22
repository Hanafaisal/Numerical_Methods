#include "integration.h"

#include <cmath>
#include <stdexcept>

namespace numerical {

namespace {

constexpr double kUniformTolerance = 1e-6;

void requireSampled(const std::vector<double>& x, const std::vector<double>& y,
                    std::size_t minPoints) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("integration: x and y must have the same length");
    }
    if (x.size() < minPoints) {
        throw std::invalid_argument("integration: not enough points");
    }
}

}  // namespace

double trapezoid(const std::vector<double>& x, const std::vector<double>& y) {
    requireSampled(x, y, 2);
    double sum = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        const double dx = x[i + 1] - x[i];
        if (dx <= 0.0) {
            throw std::invalid_argument("integration: x must be strictly increasing");
        }
        sum += 0.5 * dx * (y[i] + y[i + 1]);
    }
    return sum;
}

double simpson(const std::vector<double>& x, const std::vector<double>& y) {
    requireSampled(x, y, 3);
    const std::size_t n = x.size();
    if (n % 2 == 0) {
        throw std::invalid_argument(
            "integration: Simpson's rule needs an odd number of points (even number of intervals)");
    }

    const double h = (x.back() - x.front()) / static_cast<double>(n - 1);
    if (h <= 0.0) {
        throw std::invalid_argument("integration: x must be strictly increasing");
    }
    for (std::size_t i = 0; i + 1 < n; ++i) {
        if (std::fabs((x[i + 1] - x[i]) - h) > kUniformTolerance * std::fmax(1.0, h)) {
            throw std::invalid_argument("integration: Simpson's rule needs uniformly spaced x");
        }
    }

    double odd = 0.0, even = 0.0;
    for (std::size_t i = 1; i + 1 < n; ++i) {
        (i % 2 == 1 ? odd : even) += y[i];
    }
    return h / 3.0 * (y.front() + y.back() + 4.0 * odd + 2.0 * even);
}

}  // namespace numerical
