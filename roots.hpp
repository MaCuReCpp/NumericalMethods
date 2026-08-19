#pragma once

#include "types.hpp"

namespace numerical {

IterationResult bisection(
    const ScalarFunction& f,
    double left,
    double right,
    Tolerance tolerance = {},
    std::size_t max_iterations = 100);

IterationResult newton(
    const ScalarFunction& f,
    const ScalarFunction& derivative,
    double initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 50);

IterationResult secant(
    const ScalarFunction& f,
    double x0,
    double x1,
    Tolerance tolerance = {},
    std::size_t max_iterations = 100);

IterationResult fixed_point_iteration(
    const ScalarFunction& g,
    double initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 100);

IterationResult steffensen(
    const ScalarFunction& g,
    double initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 100);

Complex muller(
    const ScalarFunction& f,
    double x0,
    double x1,
    double x2,
    Tolerance tolerance = {},
    std::size_t max_iterations = 100);

} // namespace numerical
