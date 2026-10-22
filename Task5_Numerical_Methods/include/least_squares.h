#ifndef LEAST_SQUARES_H
#define LEAST_SQUARES_H

#include <cstddef>
#include <vector>

namespace numerical {

// Result of a polynomial least-squares fit  y ~ c0 + c1 x + c2 x^2 + ...
struct PolynomialFit {
    std::vector<double> coefficients;  // coefficients[k] multiplies x^k
    double rSquared;                   // coefficient of determination
};

// Fits a polynomial of the given degree by solving the normal equations
// (Gaussian elimination with partial pivoting).
// Throws std::invalid_argument if x and y differ in length or there are fewer
// than degree + 1 points, and std::runtime_error if the system is singular.
PolynomialFit fitPolynomial(const std::vector<double>& x,
                            const std::vector<double>& y,
                            std::size_t degree);

// Evaluates a polynomial given by its coefficients (Horner's scheme).
double evaluatePolynomial(const std::vector<double>& coefficients, double x);

}  // namespace numerical

#endif  // LEAST_SQUARES_H
