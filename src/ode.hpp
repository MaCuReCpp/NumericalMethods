#pragma once

#include "types.hpp"

#include <vector>

namespace numerical {

std::vector<OdePoint> euler_method(
    const OdeFunction& f,
    double t0,
    double y0,
    double step,
    std::size_t steps);

std::vector<OdePoint> runge_kutta4(
    const OdeFunction& f,
    double t0,
    double y0,
    double step,
    std::size_t steps);

std::vector<OdePoint> runge_kutta_fehlberg45(
    const OdeFunction& f,
    double t0,
    double y0,
    double t_end,
    double initial_step,
    double tolerance = 1e-6,
    double min_step = 1e-8,
    double max_step = 0.1);

std::vector<OdePoint> adams_bashforth_moulton4(
    const OdeFunction& f,
    double t0,
    double y0,
    double step,
    std::size_t steps);

Grid1D linear_shooting_bvp(
    const ScalarFunction& p,
    const ScalarFunction& q,
    const ScalarFunction& r,
    double a,
    double b,
    double alpha,
    double beta,
    std::size_t steps);

} // namespace numerical
