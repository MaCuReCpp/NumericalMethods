#pragma once

#include "types.hpp"

namespace numerical {

double finite_difference(
    const ScalarFunction& f,
    double x,
    double step = 1e-5,
    DifferenceMethod method = DifferenceMethod::Central);

double richardson_derivative(
    const ScalarFunction& f,
    double x,
    double initial_step,
    std::size_t levels);

double trapezoidal_rule(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t subintervals);

double simpson_rule(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t subintervals);

double romberg_integration(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t levels);

double adaptive_simpson(
    const ScalarFunction& f,
    double a,
    double b,
    double tolerance = 1e-8,
    std::size_t max_depth = 20);

double gaussian_quadrature(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t points);

} // namespace numerical
