#ifndef MOVING_AVERAGE_H
#define MOVING_AVERAGE_H

#include <cstddef>
#include <vector>

namespace numerical {

// Centered moving-average (low-pass) filter.
//
// out[i] is the mean of x[i-h .. i+h] with h = window / 2. Near the ends the
// window is truncated to the samples that exist, so the output has the same
// length as the input and no phase lag.
//
// `window` must be odd and >= 1, otherwise std::invalid_argument is thrown.
std::vector<double> movingAverage(const std::vector<double>& x, std::size_t window);

}  // namespace numerical

#endif  // MOVING_AVERAGE_H
