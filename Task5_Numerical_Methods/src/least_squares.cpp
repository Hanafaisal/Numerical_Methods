#include "least_squares.h"

#include <cmath>
#include <stdexcept>
#include <utility>

namespace numerical {

namespace {

// Solves A c = b in place (A is m x m) with partial pivoting.
std::vector<double> solveLinearSystem(std::vector<std::vector<double>> A, std::vector<double> b) {
    const std::size_t m = b.size();
    for (std::size_t col = 0; col < m; ++col) {
        std::size_t pivot = col;
        for (std::size_t r = col + 1; r < m; ++r) {
            if (std::fabs(A[r][col]) > std::fabs(A[pivot][col])) {
                pivot = r;
            }
        }
        if (std::fabs(A[pivot][col]) < 1e-12) {
            throw std::runtime_error("least_squares: normal equations are singular");
        }
        std::swap(A[col], A[pivot]);
        std::swap(b[col], b[pivot]);

        for (std::size_t r = col + 1; r < m; ++r) {
            const double f = A[r][col] / A[col][col];
            for (std::size_t c = col; c < m; ++c) {
                A[r][c] -= f * A[col][c];
            }
            b[r] -= f * b[col];
        }
    }

    std::vector<double> sol(m);
    for (std::size_t i = m; i-- > 0;) {
        double s = b[i];
        for (std::size_t c = i + 1; c < m; ++c) {
            s -= A[i][c] * sol[c];
        }
        sol[i] = s / A[i][i];
    }
    return sol;
}

}  // namespace

double evaluatePolynomial(const std::vector<double>& coefficients, double x) {
    double result = 0.0;
    for (std::size_t k = coefficients.size(); k-- > 0;) {
        result = result * x + coefficients[k];
    }
    return result;
}

PolynomialFit fitPolynomial(const std::vector<double>& x, const std::vector<double>& y,
                            std::size_t degree) {
    if (x.size() != y.size()) {
        throw std::invalid_argument("least_squares: x and y must have the same length");
    }
    if (x.size() < degree + 1) {
        throw std::invalid_argument("least_squares: need at least degree + 1 points");
    }

    const std::size_t m = degree + 1;

    // Power sums: s[k] = sum x^k (k = 0 .. 2*degree), t[k] = sum y x^k.
    std::vector<double> s(2 * degree + 1, 0.0);
    std::vector<double> t(m, 0.0);
    for (std::size_t i = 0; i < x.size(); ++i) {
        double p = 1.0;
        for (std::size_t k = 0; k <= 2 * degree; ++k) {
            s[k] += p;
            if (k < m) {
                t[k] += y[i] * p;
            }
            p *= x[i];
        }
    }

    std::vector<std::vector<double>> A(m, std::vector<double>(m));
    for (std::size_t i = 0; i < m; ++i) {
        for (std::size_t j = 0; j < m; ++j) {
            A[i][j] = s[i + j];
        }
    }

    PolynomialFit fit;
    fit.coefficients = solveLinearSystem(A, t);

    double mean = 0.0;
    for (double v : y) {
        mean += v;
    }
    mean /= static_cast<double>(y.size());

    double ssRes = 0.0, ssTot = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        const double r = y[i] - evaluatePolynomial(fit.coefficients, x[i]);
        ssRes += r * r;
        ssTot += (y[i] - mean) * (y[i] - mean);
    }
    fit.rSquared = ssTot > 0.0 ? 1.0 - ssRes / ssTot : 1.0;
    return fit;
}

}  // namespace numerical
