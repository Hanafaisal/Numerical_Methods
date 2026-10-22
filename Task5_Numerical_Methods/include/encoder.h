#ifndef ENCODER_H
#define ENCODER_H

#include <cstdint>
#include <vector>

namespace numerical {

struct EncoderConfig {
    double countsPerRevolution = 1024.0;  // counts per revolution after quadrature decoding
    unsigned counterBits = 16;            // width of the hardware counter, 2..32
};

struct EncoderSample {
    double time;                     // seconds
    std::int64_t cumulativeCounts;   // counts since the first sample, counter wrap removed
    double angleDeg;                 // degrees since the first sample
    double velocityDegPerSec;        // backward difference; NaN for the first sample
    double velocityRpm;              // same velocity in revolutions per minute
};

// Converts raw wrapping counter readings into angle and angular velocity.
//
// The counter may overflow/underflow between samples; each step is interpreted
// as the shortest signed distance modulo 2^counterBits, which is correct as long
// as the shaft moves less than half the counter range between two samples.
//
// Throws std::invalid_argument for mismatched/empty inputs, a raw count that does
// not fit in counterBits, an invalid config, or non-increasing time.
std::vector<EncoderSample> decodeEncoder(const std::vector<double>& time,
                                         const std::vector<std::uint32_t>& rawCounts,
                                         const EncoderConfig& config);

}  // namespace numerical

#endif  // ENCODER_H
