#pragma once

#include "types.hpp"

namespace numerical {

Vector least_squares_polynomial(
    const Vector& xs,
    const Vector& ys,
    std::size_t degree);

double evaluate_polynomial(
    const Vector& coefficients,
    double x);

Vector chebyshev_nodes(double a, double b, std::size_t count);
double chebyshev_polynomial(int degree, double x);
ComplexVector fft(const ComplexVector& values, bool inverse = false);

} // namespace numerical
