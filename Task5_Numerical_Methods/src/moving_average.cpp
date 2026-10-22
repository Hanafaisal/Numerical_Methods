#include "moving_average.h"

#include <algorithm>
#include <stdexcept>

namespace numerical {

std::vector<double> movingAverage(const std::vector<double>& x, std::size_t window) {
    if (window == 0 || window % 2 == 0) {
        throw std::invalid_argument("moving_average: window must be odd and >= 1");
    }

    const std::size_t n = x.size();
    std::vector<double> out(n);
    if (n == 0) {
        return out;
    }

    // Prefix sums make every window mean O(1).
    std::vector<double> prefix(n + 1, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        prefix[i + 1] = prefix[i] + x[i];
    }

    const std::size_t half = window / 2;
    for (std::size_t i = 0; i < n; ++i) {
        const std::size_t lo = (i >= half) ? i - half : 0;
        const std::size_t hi = std::min(i + half, n - 1);
        out[i] = (prefix[hi + 1] - prefix[lo]) / static_cast<double>(hi - lo + 1);
    }
    return out;
}

}  // namespace numerical
