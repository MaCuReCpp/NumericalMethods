#include "nonlinear.hpp"
#include "roots.hpp"
#include "internal.hpp"
#include "linear_algebra.hpp"

#include <cmath>
#include <limits>

namespace numerical {
using detail::l2_distance;
using detail::subtract;

IterationResult fixed_point_iteration(
    const ScalarFunction& g,
    double initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    double previous = initial_guess;
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const double current = g(previous);
        const double scale = std::max(1.0, std::abs(current));
        const double error = std::abs(current - previous);
        if (error <= tolerance.absolute + tolerance.relative * scale) {
            return {current, error, iteration, true};
        }
        previous = current;
    }
    return {previous, std::numeric_limits<double>::quiet_NaN(), max_iterations, false};
}

IterationResult steffensen(
    const ScalarFunction& g,
    double initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    double x = initial_guess;
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const double x1 = g(x);
        const double x2 = g(x1);
        const double denominator = x2 - 2.0 * x1 + x;
        if (std::abs(denominator) < detail::kTiny) {
            throw std::runtime_error("Steffensen method encountered a near-zero denominator");
        }
        const double accelerated = x - (x1 - x) * (x1 - x) / denominator;
        const double error = std::abs(accelerated - x);
        const double scale = std::max(1.0, std::abs(accelerated));
        if (error <= tolerance.absolute + tolerance.relative * scale) {
            return {accelerated, error, iteration, true};
        }
        x = accelerated;
    }
    return {x, std::numeric_limits<double>::quiet_NaN(), max_iterations, false};
}

Complex muller(
    const ScalarFunction& f,
    double x0,
    double x1,
    double x2,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    Complex z0 = x0;
    Complex z1 = x1;
    Complex z2 = x2;
    auto fc = [&](Complex z) { return Complex(f(z.real()), 0.0); };

    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        const Complex h0 = z1 - z0;
        const Complex h1 = z2 - z1;
        const Complex delta0 = (fc(z1) - fc(z0)) / h0;
        const Complex delta1 = (fc(z2) - fc(z1)) / h1;
        const Complex a = (delta1 - delta0) / (h1 + h0);
        const Complex b = a * h1 + delta1;
        const Complex c = fc(z2);
        const Complex discriminant = std::sqrt(b * b - 4.0 * a * c);
        const Complex denominator = (std::abs(b + discriminant) > std::abs(b - discriminant))
            ? b + discriminant
            : b - discriminant;
        if (std::abs(denominator) < detail::kTiny) {
            throw std::runtime_error("Muller method encountered a near-zero denominator");
        }

        const Complex step = -2.0 * c / denominator;
        const Complex next = z2 + step;
        if (std::abs(step) <= tolerance.absolute + tolerance.relative * std::max(1.0, std::abs(next))) {
            return next;
        }

        z0 = z1;
        z1 = z2;
        z2 = next;
    }

    return z2;
}

LinearSolveResult nonlinear_newton(
    const VectorFunction& f,
    const JacobianFunction& jacobian,
    Vector initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    Vector x = initial_guess;
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const Vector fx = f(x);
        const double residual = l2_norm(fx);
        if (residual <= tolerance.absolute) {
            return {x, residual, iteration - 1, true};
        }
        Vector rhs(fx.size(), 0.0);
        for (std::size_t i = 0; i < fx.size(); ++i) {
            rhs[i] = -fx[i];
        }
        const Vector delta = gaussian_elimination(jacobian(x), rhs).solution;
        Vector next = x;
        for (std::size_t i = 0; i < x.size(); ++i) {
            next[i] += delta[i];
        }
        if (l2_distance(next, x) <= tolerance.absolute + tolerance.relative * std::max(1.0, l2_norm(next))) {
            return {next, l2_norm(f(next)), iteration, true};
        }
        x = next;
    }
    return {x, l2_norm(f(x)), max_iterations, false};
}

} // namespace numerical
