# Task 5 – Numerical Methods Pipeline

A small C++17 project that runs four numerical-methods stages over CSV data and
writes the results back to CSV:

| Stage | Module | Input | Output |
|-------|--------|-------|--------|
| 1. Error analysis + least-squares calibration | `error_analysis`, `least_squares` | `reference_data.csv`, `measured_data.csv` | `results/error_results.csv` |
| 2. Moving-average filtering | `moving_average` | `noisy_data.csv` | `results/filtered_data.csv` |
| 3. Numerical integration | `integration` | reference, measured, noisy and filtered series | `results/integration_results.csv` |
| 4. Encoder decoding + least-squares speed fit | `encoder`, `least_squares` | `encoder_data.csv` | `results/encoder_results.csv` |

## Project layout

```
Task5_Numerical_Methods/
├── main.cpp                 pipeline driver (CSV I/O, runs the 4 stages)
├── src/                     module implementations
├── include/                 module headers (namespace numerical)
├── data/                    input CSVs
├── results/                 output CSVs (regenerated on every run)
└── README.md
```

## Build and run

Requires a C++17 compiler (tested with g++ 13).

```bash
cd Task5_Numerical_Methods
g++ -std=c++17 -O2 -Wall -Wextra -Iinclude main.cpp src/*.cpp -o numerical_pipeline
./numerical_pipeline            # uses the current directory as project root
./numerical_pipeline /path/to/Task5_Numerical_Methods   # or give the root explicitly
```

The program reads `data/*.csv`, creates `results/` if needed, prints a summary
to the console and overwrites the four result files. It exits with status 1 and
an `error:` message on bad input (missing file, malformed row, mismatched
lengths, …).

## Modules

**`error_analysis`** – `signedErrors`, `absoluteErrors`, `relativeErrorsPercent`,
`computeErrorStats` (bias, MAE, MSE, RMSE, max |error|, mean relative error).
Relative error is undefined where the reference is 0; those samples are left
blank in the CSV and excluded from the mean.

**`moving_average`** – `movingAverage(x, window)`: centered, odd window, O(n)
using prefix sums. Windows are truncated at both ends so the output has the same
length as the input and no phase lag.

**`integration`** – `trapezoid(x, y)` (any strictly increasing x) and
`simpson(x, y)` (uniform spacing, odd number of points – it throws otherwise).

**`least_squares`** – `fitPolynomial(x, y, degree)` solves the normal equations
with partial pivoting and returns the coefficients and R²;
`evaluatePolynomial` uses Horner's scheme.

**`encoder`** – `decodeEncoder(time, rawCounts, config)` turns raw counter
readings into cumulative counts, angle (deg), and angular velocity (deg/s and
rpm). Counter overflow/underflow is handled by taking the shortest signed step
modulo 2^bits, which is valid while the shaft moves less than half the counter
range between samples. Velocity is a backward difference, so the first sample
has none (blank in the CSV).

## Data files

All inputs are two-column CSVs with a header row.

| File | Columns | Content (sample data) |
|------|---------|-----------------------|
| `reference_data.csv` | `time,value` | true signal `3·sin(0.5t) + 0.5t`, 101 samples, dt = 0.1 s |
| `measured_data.csv` | `time,value` | reference with +3 % gain, +0.20 offset and Gaussian noise (σ = 0.08) |
| `noisy_data.csv` | `time,value` | `2·sin(1.5t) + 1` plus Gaussian noise (σ = 0.4), 201 samples, dt = 0.05 s |
| `encoder_data.csv` | `time,count` | 16-bit counter, 1024 counts/rev, starts at 65400 so it wraps; speed `120 + 30·sin(0.8t)` deg/s |

The sample data is synthetic and generated with a fixed seed. **If you replace it
with real data**, note that `main.cpp` contains the known underlying functions
(`referenceAntiderivative`, `cleanNoisySignal`, `noisyAntiderivative`) used to
report "exact" integrals and error against the clean signal. Update or remove
those for data that has no known ground truth. The encoder resolution and counter
width are the constants `kCountsPerRev` and `kCounterBits` at the top of `main.cpp`,
and the filter window is `kFilterWindow`.

## Result files

- **`error_results.csv`** – `time, reference, measured, error, abs_error,
  rel_error_percent, corrected, corrected_error`. `corrected` is the measurement
  after inverting the least-squares calibration `measured ≈ offset + gain·reference`.
- **`filtered_data.csv`** – `time, noisy, filtered`.
- **`integration_results.csv`** – `dataset, method, integral, exact_integral,
  abs_error, rel_error_percent` for trapezoid and Simpson on each of the
  reference, measured, noisy and filtered series.
- **`encoder_results.csv`** – `time, raw_count, cumulative_counts, angle_deg,
  velocity_deg_s, velocity_rpm, fitted_angle_deg` (the last column is the
  least-squares line through angle vs time; its slope is the average speed).

## Sample run

```
[1] Error analysis (101 samples)
  raw measurement       bias   0.2932  MAE  0.2932  RMSE  0.3036  max|e|  0.5088  mean rel% 12.8479
  after calibration     bias  -0.0000  MAE  0.0552  RMSE  0.0689  max|e|  0.2046  mean rel%  2.3924
  least-squares fit: measured = 0.2151 + 1.0268 * reference   (R^2 = 0.99717)

[2] Moving-average filter (window = 7, 201 samples)
  noisy vs clean        RMSE  0.4218
  filtered vs clean     RMSE  0.1616

[3] Numerical integration
  reference trapezoid    29.29713  (exact 29.29803)
  reference simpson      29.29803  (exact 29.29803)
  noisy     trapezoid    12.73680  (exact 12.34625)
  filtered  trapezoid    12.74089  (exact 12.34625)

[4] Encoder
  counter wrap-arounds handled: 1
  total rotation: 1242.77 deg (3.452 rev)
  average speed (least-squares slope): 118.114 deg/s = 19.686 rpm
```

Things worth noticing in these numbers:

- Calibration removes the systematic part of the sensor error (bias → 0, RMSE
  0.30 → 0.07); what remains is the random noise.
- The filter cuts the RMSE against the clean signal by about 60 %, but the
  integral of the filtered signal is *not* more accurate than the raw one.
  Integration already averages out zero-mean noise, and the filter slightly
  flattens the peaks. The remaining integral error comes from noise that
  happens not to cancel in this particular sample.
- Simpson beats the trapezoid rule on the smooth reference signal, but on noisy
  data the two are limited by the noise, not the method.
- The measured series integrates ~10 % too high because of its gain and offset
  error, exactly what the calibration in stage 1 estimates.
