#pragma once

#include "types.hpp"

namespace numerical {

LinearSolveResult nonlinear_newton(
    const VectorFunction& f,
    const JacobianFunction& jacobian,
    Vector initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 50);

} // namespace numerical
