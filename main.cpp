// Task 5 - Numerical Methods pipeline
//
//   data/*.csv  ->  error analysis (+ least-squares calibration)
//               ->  moving-average filtering
//               ->  numerical integration (trapezoid / Simpson)
//               ->  encoder decoding (+ least-squares speed fit)
//               ->  results/*.csv
//
// Usage: numerical_pipeline [project_root]     (default: current directory)

#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "encoder.h"
#include "error_analysis.h"
#include "integration.h"
#include "least_squares.h"
#include "moving_average.h"

namespace fs = std::filesystem;

namespace {

// ---- configuration -------------------------------------------------------
constexpr std::size_t kFilterWindow = 7;       // moving-average window (odd)
constexpr double kCountsPerRev = 1024.0;       // encoder counts per revolution
constexpr unsigned kCounterBits = 16;          // encoder counter width
constexpr double kTimeTolerance = 1e-6;        // for matching time stamps

// ---- known "ground truth" of the sample data (see README) ----------------
// reference_data.csv / measured_data.csv are based on y(t) = 3 sin(0.5 t) + 0.5 t
double referenceAntiderivative(double x) { return -6.0 * std::cos(0.5 * x) + 0.25 * x * x; }
// noisy_data.csv is based on y(t) = 2 sin(1.5 t) + 1
double cleanNoisySignal(double t) { return 2.0 * std::sin(1.5 * t) + 1.0; }
double noisyAntiderivative(double x) { return -(2.0 / 1.5) * std::cos(1.5 * x) + x; }

// ---- CSV helpers ---------------------------------------------------------
struct Series {
    std::vector<double> t;
    std::vector<double> y;
};

Series readSeries(const fs::path& path) {
    std::ifstream in(path);
    if (!in) {
        throw std::runtime_error("cannot open " + path.string());
    }

    Series s;
    std::string line;
    std::size_t lineNo = 0;
    while (std::getline(in, line)) {
        ++lineNo;
        if (lineNo == 1 || line.empty() || line == "\r") {
            continue;  // header / blank line
        }
        std::stringstream ss(line);
        std::string a, b;
        if (!std::getline(ss, a, ',') || !std::getline(ss, b, ',')) {
            throw std::runtime_error(path.string() + ": expected 2 columns on line " +
                                     std::to_string(lineNo));
        }
        try {
            s.t.push_back(std::stod(a));
            s.y.push_back(std::stod(b));
        } catch (const std::exception&) {
            throw std::runtime_error(path.string() + ": invalid number on line " +
                                     std::to_string(lineNo));
        }
    }
    if (s.t.empty()) {
        throw std::runtime_error(path.string() + ": no data rows");
    }
    return s;
}

// NaN is written as an empty field.
std::string num(double v, int precision = 6) {
    if (std::isnan(v)) {
        return "";
    }
    std::ostringstream o;
    o << std::fixed << std::setprecision(precision) << v;
    return o.str();
}

using Row = std::vector<std::string>;

void writeCsv(const fs::path& path, const std::string& header, const std::vector<Row>& rows) {
    std::ofstream out(path);
    if (!out) {
        throw std::runtime_error("cannot write " + path.string());
    }
    out << header << '\n';
    for (const Row& r : rows) {
        for (std::size_t i = 0; i < r.size(); ++i) {
            out << (i ? "," : "") << r[i];
        }
        out << '\n';
    }
}

void requireSameTimes(const Series& a, const Series& b, const std::string& what) {
    if (a.t.size() != b.t.size()) {
        throw std::runtime_error(what + ": the two files have different numbers of samples");
    }
    for (std::size_t i = 0; i < a.t.size(); ++i) {
        if (std::fabs(a.t[i] - b.t[i]) > kTimeTolerance) {
            throw std::runtime_error(what + ": time stamps differ at sample " + std::to_string(i));
        }
    }
}

void printStats(const std::string& label, const numerical::ErrorStats& s) {
    std::cout << "  " << std::left << std::setw(22) << label << std::right << std::fixed
              << std::setprecision(4) << "bias " << std::setw(8) << s.bias << "  MAE "
              << std::setw(7) << s.mae << "  RMSE " << std::setw(7) << s.rmse << "  max|e| "
              << std::setw(7) << s.maxAbsError << "  mean rel% " << std::setw(7)
              << s.meanRelPercent << '\n';
}

// ---- stage 1: error analysis ---------------------------------------------
void runErrorAnalysis(const fs::path& results, const Series& ref, const Series& meas) {
    requireSameTimes(ref, meas, "reference/measured");

    const auto err = numerical::signedErrors(ref.y, meas.y);
    const auto absErr = numerical::absoluteErrors(ref.y, meas.y);
    const auto relErr = numerical::relativeErrorsPercent(ref.y, meas.y);
    const auto before = numerical::computeErrorStats(ref.y, meas.y);

    // Calibrate the sensor with least squares: measured ~ offset + gain * reference,
    // then invert the fit to correct the measurements.
    const auto fit = numerical::fitPolynomial(ref.y, meas.y, 1);
    const double offset = fit.coefficients[0];
    const double gain = fit.coefficients[1];
    if (std::fabs(gain) < 1e-12) {
        throw std::runtime_error("calibration gain is ~0, cannot correct measurements");
    }
    std::vector<double> corrected(meas.y.size());
    for (std::size_t i = 0; i < corrected.size(); ++i) {
        corrected[i] = (meas.y[i] - offset) / gain;
    }
    const auto after = numerical::computeErrorStats(ref.y, corrected);

    std::vector<Row> rows;
    for (std::size_t i = 0; i < ref.t.size(); ++i) {
        rows.push_back({num(ref.t[i], 4), num(ref.y[i]), num(meas.y[i]), num(err[i]),
                        num(absErr[i]), num(relErr[i]), num(corrected[i]),
                        num(corrected[i] - ref.y[i])});
    }
    writeCsv(results / "error_results.csv",
             "time,reference,measured,error,abs_error,rel_error_percent,corrected,corrected_error",
             rows);

    std::cout << "[1] Error analysis (" << ref.t.size() << " samples)\n";
    printStats("raw measurement", before);
    printStats("after calibration", after);
    std::cout << "  least-squares fit: measured = " << num(offset, 4) << " + " << num(gain, 4)
              << " * reference   (R^2 = " << num(fit.rSquared, 5) << ")\n"
              << "  relative error skipped for " << ref.t.size() - before.relativeSamples
              << " sample(s) where reference = 0\n";
}

// ---- stage 2: moving-average filter --------------------------------------
std::vector<double> runFiltering(const fs::path& results, const Series& noisy) {
    const std::vector<double> filtered = numerical::movingAverage(noisy.y, kFilterWindow);

    std::vector<Row> rows;
    for (std::size_t i = 0; i < noisy.t.size(); ++i) {
        rows.push_back({num(noisy.t[i], 4), num(noisy.y[i]), num(filtered[i])});
    }
    writeCsv(results / "filtered_data.csv", "time,noisy,filtered", rows);

    // The sample signal is synthetic, so we can measure how close each version is to the truth.
    std::vector<double> clean(noisy.t.size());
    for (std::size_t i = 0; i < clean.size(); ++i) {
        clean[i] = cleanNoisySignal(noisy.t[i]);
    }
    std::cout << "[2] Moving-average filter (window = " << kFilterWindow << ", "
              << noisy.t.size() << " samples)\n";
    printStats("noisy vs clean", numerical::computeErrorStats(clean, noisy.y));
    printStats("filtered vs clean", numerical::computeErrorStats(clean, filtered));
    return filtered;
}

// ---- stage 3: numerical integration --------------------------------------
struct Dataset {
    std::string name;
    std::vector<double> t;
    std::vector<double> y;
    double exact;  // exact integral of the underlying true signal
};

void runIntegration(const fs::path& results, const Series& ref, const Series& meas,
                    const Series& noisy, const std::vector<double>& filtered) {
    const double refLo = ref.t.front(), refHi = ref.t.back();
    const double noisyLo = noisy.t.front(), noisyHi = noisy.t.back();
    const double exactRef = referenceAntiderivative(refHi) - referenceAntiderivative(refLo);
    const double exactNoisy = noisyAntiderivative(noisyHi) - noisyAntiderivative(noisyLo);

    const std::vector<Dataset> datasets = {
        {"reference", ref.t, ref.y, exactRef},
        {"measured", meas.t, meas.y, exactRef},
        {"noisy", noisy.t, noisy.y, exactNoisy},
        {"filtered", noisy.t, filtered, exactNoisy},
    };

    std::vector<Row> rows;
    std::cout << "[3] Numerical integration\n";
    for (const Dataset& d : datasets) {
        const double results2[2] = {numerical::trapezoid(d.t, d.y), numerical::simpson(d.t, d.y)};
        const char* methods[2] = {"trapezoid", "simpson"};
        for (int m = 0; m < 2; ++m) {
            const double absError = std::fabs(results2[m] - d.exact);
            rows.push_back({d.name, methods[m], num(results2[m]), num(d.exact), num(absError),
                            num(absError / std::fabs(d.exact) * 100.0)});
            std::cout << "  " << std::left << std::setw(10) << d.name << std::setw(11)
                      << methods[m] << std::right << std::fixed << std::setprecision(5)
                      << std::setw(10) << results2[m] << "  (exact " << d.exact << ", |err| "
                      << absError << ")\n";
        }
    }
    writeCsv(results / "integration_results.csv",
             "dataset,method,integral,exact_integral,abs_error,rel_error_percent", rows);
}

// ---- stage 4: encoder ----------------------------------------------------
void runEncoder(const fs::path& results, const Series& enc) {
    std::vector<std::uint32_t> raw(enc.y.size());
    for (std::size_t i = 0; i < raw.size(); ++i) {
        if (enc.y[i] < 0.0) {
            throw std::runtime_error("encoder_data.csv: negative count at sample " +
                                     std::to_string(i));
        }
        raw[i] = static_cast<std::uint32_t>(std::llround(enc.y[i]));
    }

    numerical::EncoderConfig cfg;
    cfg.countsPerRevolution = kCountsPerRev;
    cfg.counterBits = kCounterBits;
    const auto samples = numerical::decodeEncoder(enc.t, raw, cfg);

    std::vector<double> angle(samples.size());
    for (std::size_t i = 0; i < samples.size(); ++i) {
        angle[i] = samples[i].angleDeg;
    }
    // Average angular velocity = slope of the least-squares line angle(t).
    const auto fit = numerical::fitPolynomial(enc.t, angle, 1);

    std::vector<Row> rows;
    std::size_t wraps = 0;
    for (std::size_t i = 0; i < samples.size(); ++i) {
        const auto& s = samples[i];
        if (i > 0 && raw[i] < raw[i - 1] && s.cumulativeCounts > samples[i - 1].cumulativeCounts) {
            ++wraps;
        }
        rows.push_back({num(s.time, 4), std::to_string(raw[i]), std::to_string(s.cumulativeCounts),
                        num(s.angleDeg), num(s.velocityDegPerSec), num(s.velocityRpm),
                        num(numerical::evaluatePolynomial(fit.coefficients, s.time))});
    }
    writeCsv(results / "encoder_results.csv",
             "time,raw_count,cumulative_counts,angle_deg,velocity_deg_s,velocity_rpm,"
             "fitted_angle_deg",
             rows);

    std::cout << "[4] Encoder (" << samples.size() << " samples, " << num(kCountsPerRev, 0)
              << " counts/rev, " << kCounterBits << "-bit counter)\n"
              << "  counter wrap-arounds handled: " << wraps << '\n'
              << "  total rotation: " << num(samples.back().angleDeg, 2) << " deg ("
              << num(samples.back().angleDeg / 360.0, 3) << " rev)\n"
              << "  average speed (least-squares slope): " << num(fit.coefficients[1], 3)
              << " deg/s = " << num(fit.coefficients[1] / 6.0, 3)
              << " rpm   (R^2 = " << num(fit.rSquared, 5) << ")\n";
}

}  // namespace

int main(int argc, char** argv) {
    try {
        const fs::path root = argc > 1 ? fs::path(argv[1]) : fs::current_path();
        const fs::path data = root / "data";
        const fs::path results = root / "results";
        fs::create_directories(results);

        const Series ref = readSeries(data / "reference_data.csv");
        const Series meas = readSeries(data / "measured_data.csv");
        const Series noisy = readSeries(data / "noisy_data.csv");
        const Series enc = readSeries(data / "encoder_data.csv");

        std::cout << "Numerical methods pipeline  (root: " << root.string() << ")\n\n";
        runErrorAnalysis(results, ref, meas);
        std::cout << '\n';
        const std::vector<double> filtered = runFiltering(results, noisy);
        std::cout << '\n';
        runIntegration(results, ref, meas, noisy, filtered);
        std::cout << '\n';
        runEncoder(results, enc);
        std::cout << "\nDone. Results written to " << results.string() << '\n';
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << '\n';
        return 1;
    }
}
