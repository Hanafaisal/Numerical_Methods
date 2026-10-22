#include "error_analysis.h"

#include <cmath>
#include <limits>
#include <stdexcept>

namespace numerical {

namespace {

constexpr double kZeroReference = 1e-12;

void requireValid(const std::vector<double>& reference, const std::vector<double>& measured) {
    if (reference.empty() || reference.size() != measured.size()) {
        throw std::invalid_argument(
            "error_analysis: inputs must be non-empty and have the same length");
    }
}

}  // namespace

std::vector<double> signedErrors(const std::vector<double>& reference,
                                 const std::vector<double>& measured) {
    requireValid(reference, measured);
    std::vector<double> out(reference.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        out[i] = measured[i] - reference[i];
    }
    return out;
}

std::vector<double> absoluteErrors(const std::vector<double>& reference,
                                   const std::vector<double>& measured) {
    std::vector<double> out = signedErrors(reference, measured);
    for (double& e : out) {
        e = std::fabs(e);
    }
    return out;
}

std::vector<double> relativeErrorsPercent(const std::vector<double>& reference,
                                          const std::vector<double>& measured) {
    requireValid(reference, measured);
    std::vector<double> out(reference.size());
    for (std::size_t i = 0; i < out.size(); ++i) {
        if (std::fabs(reference[i]) < kZeroReference) {
            out[i] = std::numeric_limits<double>::quiet_NaN();
        } else {
            out[i] = std::fabs(measured[i] - reference[i]) / std::fabs(reference[i]) * 100.0;
        }
    }
    return out;
}

ErrorStats computeErrorStats(const std::vector<double>& reference,
                             const std::vector<double>& measured) {
    const std::vector<double> err = signedErrors(reference, measured);
    const std::vector<double> rel = relativeErrorsPercent(reference, measured);

    double sum = 0.0, sumAbs = 0.0, sumSq = 0.0, maxAbs = 0.0, sumRel = 0.0;
    std::size_t relCount = 0;
    for (std::size_t i = 0; i < err.size(); ++i) {
        const double a = std::fabs(err[i]);
        sum += err[i];
        sumAbs += a;
        sumSq += err[i] * err[i];
        if (a > maxAbs) {
            maxAbs = a;
        }
        if (!std::isnan(rel[i])) {
            sumRel += rel[i];
            ++relCount;
        }
    }

    const double n = static_cast<double>(err.size());
    ErrorStats s;
    s.bias = sum / n;
    s.mae = sumAbs / n;
    s.mse = sumSq / n;
    s.rmse = std::sqrt(s.mse);
    s.maxAbsError = maxAbs;
    s.relativeSamples = relCount;
    s.meanRelPercent =
        relCount > 0 ? sumRel / static_cast<double>(relCount)
                     : std::numeric_limits<double>::quiet_NaN();
    return s;
}

}  // namespace numerical
