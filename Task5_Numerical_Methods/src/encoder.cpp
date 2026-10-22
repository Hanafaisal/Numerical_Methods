#include "encoder.h"

#include <limits>
#include <stdexcept>

namespace numerical {

std::vector<EncoderSample> decodeEncoder(const std::vector<double>& time,
                                         const std::vector<std::uint32_t>& rawCounts,
                                         const EncoderConfig& config) {
    if (time.empty() || time.size() != rawCounts.size()) {
        throw std::invalid_argument("encoder: time and counts must be non-empty and equal length");
    }
    if (config.countsPerRevolution <= 0.0 || config.counterBits < 2 || config.counterBits > 32) {
        throw std::invalid_argument("encoder: invalid configuration");
    }

    const std::int64_t modulus = std::int64_t{1} << config.counterBits;
    const std::int64_t half = modulus / 2;
    const double degPerCount = 360.0 / config.countsPerRevolution;

    for (std::uint32_t c : rawCounts) {
        if (static_cast<std::int64_t>(c) >= modulus) {
            throw std::invalid_argument("encoder: raw count exceeds the counter width");
        }
    }

    std::vector<EncoderSample> out(time.size());
    out[0] = {time[0], 0, 0.0, std::numeric_limits<double>::quiet_NaN(),
              std::numeric_limits<double>::quiet_NaN()};

    for (std::size_t i = 1; i < time.size(); ++i) {
        const double dt = time[i] - time[i - 1];
        if (dt <= 0.0) {
            throw std::invalid_argument("encoder: time must be strictly increasing");
        }

        // Shortest signed step modulo the counter range (handles wrap-around).
        std::int64_t delta = static_cast<std::int64_t>(rawCounts[i]) -
                             static_cast<std::int64_t>(rawCounts[i - 1]);
        delta = ((delta % modulus) + modulus) % modulus;
        if (delta > half) {
            delta -= modulus;
        }

        EncoderSample& s = out[i];
        s.time = time[i];
        s.cumulativeCounts = out[i - 1].cumulativeCounts + delta;
        s.angleDeg = static_cast<double>(s.cumulativeCounts) * degPerCount;
        s.velocityDegPerSec = static_cast<double>(delta) * degPerCount / dt;
        s.velocityRpm = s.velocityDegPerSec / 6.0;  // 360 deg/rev, 60 s/min
    }
    return out;
}

}  // namespace numerical
