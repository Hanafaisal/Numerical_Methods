#ifndef ERROR_ANALYSIS_H
#define ERROR_ANALYSIS_H

#include <cstddef>
#include <vector>

namespace numerical {

// Summary statistics of (measured - reference).
struct ErrorStats {
    double bias;             // mean signed error
    double mae;              // mean absolute error
    double mse;              // mean squared error
    double rmse;             // root mean squared error
    double maxAbsError;      // largest |error|
    double meanRelPercent;   // mean relative error in %, NaN if no sample qualifies
    std::size_t relativeSamples;  // samples used for the relative error (|reference| != 0)
};

// All functions throw std::invalid_argument if the inputs are empty or differ in length.

// measured[i] - reference[i]
std::vector<double> signedErrors(const std::vector<double>& reference,
                                 const std::vector<double>& measured);

// |measured[i] - reference[i]|
std::vector<double> absoluteErrors(const std::vector<double>& reference,
                                   const std::vector<double>& measured);

// |measured[i] - reference[i]| / |reference[i]| * 100.
// Undefined (quiet NaN) where |reference[i]| is ~0.
std::vector<double> relativeErrorsPercent(const std::vector<double>& reference,
                                          const std::vector<double>& measured);

ErrorStats computeErrorStats(const std::vector<double>& reference,
                             const std::vector<double>& measured);

}  // namespace numerical

#endif  // ERROR_ANALYSIS_H
