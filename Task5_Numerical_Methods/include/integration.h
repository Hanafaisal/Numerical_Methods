#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <vector>

namespace numerical {

// Composite trapezoidal rule over sampled data (y vs x).
// x must be strictly increasing; spacing may be non-uniform.
// Throws std::invalid_argument for mismatched sizes, fewer than 2 points,
// or non-increasing x.
double trapezoid(const std::vector<double>& x, const std::vector<double>& y);

// Composite Simpson's 1/3 rule over sampled data.
// Requires uniformly spaced x and an odd number of points (even number of
// intervals). Throws std::invalid_argument otherwise.
double simpson(const std::vector<double>& x, const std::vector<double>& y);

}  // namespace numerical

#endif  // INTEGRATION_H
