#include "roots.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace numerical {
namespace {

constexpr double kPivotEpsilon = 1e-14;

bool close_enough(double current, double previous, double residual, Tolerance tolerance)
{
    const double scale = std::max(1.0, std::abs(current));
    const double step_error = std::abs(current - previous);
    return std::abs(residual) <= tolerance.absolute ||
           step_error <= tolerance.absolute + tolerance.relative * scale;
}

} // namespace

IterationResult bisection(
    const ScalarFunction& f,
    double left,
    double right,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    if (left >= right) {
        throw std::invalid_argument("left bound must be smaller than right bound");
    }

    double f_left = f(left);
    double f_right = f(right);
    if (f_left == 0.0) {
        return {left, 0.0, 0, true};
    }
    if (f_right == 0.0) {
        return {right, 0.0, 0, true};
    }
    if (f_left * f_right > 0.0) {
        throw std::invalid_argument("bisection requires a sign change over the interval");
    }

    double mid = left;
    double previous = mid;
    double f_mid = f_left;

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        previous = mid;
        mid = 0.5 * (left + right);
        f_mid = f(mid);

        if (close_enough(mid, previous, f_mid, tolerance) ||
            0.5 * (right - left) <= tolerance.absolute) {
            return {mid, f_mid, iteration, true};
        }

        if (f_left * f_mid < 0.0) {
            right = mid;
            f_right = f_mid;
        } else {
            left = mid;
            f_left = f_mid;
        }
    }

    (void)f_right;
    return {mid, f_mid, max_iterations, false};
}

IterationResult newton(
    const ScalarFunction& f,
    const ScalarFunction& derivative,
    double initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    double x = initial_guess;
    double residual = f(x);

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const double slope = derivative(x);
        if (std::abs(slope) < kPivotEpsilon) {
            throw std::runtime_error("Newton method encountered a near-zero derivative");
        }

        const double previous = x;
        x = x - residual / slope;
        residual = f(x);

        if (close_enough(x, previous, residual, tolerance)) {
            return {x, residual, iteration, true};
        }
    }

    return {x, residual, max_iterations, false};
}

IterationResult secant(
    const ScalarFunction& f,
    double x0,
    double x1,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    double f0 = f(x0);
    double f1 = f(x1);

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const double denominator = f1 - f0;
        if (std::abs(denominator) < kPivotEpsilon) {
            throw std::runtime_error("secant method encountered a near-zero denominator");
        }

        const double x2 = x1 - f1 * (x1 - x0) / denominator;
        const double f2 = f(x2);

        if (close_enough(x2, x1, f2, tolerance)) {
            return {x2, f2, iteration, true};
        }

        x0 = x1;
        f0 = f1;
        x1 = x2;
        f1 = f2;
    }

    return {x1, f1, max_iterations, false};
}

} // namespace numerical
