#pragma once

#include "types.hpp"

namespace numerical {

Grid1D finite_difference_bvp(
    const ScalarFunction& p,
    const ScalarFunction& q,
    const ScalarFunction& r,
    double a,
    double b,
    double alpha,
    double beta,
    std::size_t interior_points);

HeatEquationResult heat_equation_explicit(
    const ScalarFunction& initial,
    const ScalarFunction& left_boundary,
    const ScalarFunction& right_boundary,
    double x_left,
    double x_right,
    double t_final,
    std::size_t x_steps,
    std::size_t t_steps,
    double diffusivity = 1.0);

WaveEquationResult wave_equation_explicit(
    const ScalarFunction& initial_displacement,
    const ScalarFunction& initial_velocity,
    const ScalarFunction& left_boundary,
    const ScalarFunction& right_boundary,
    double x_left,
    double x_right,
    double t_final,
    std::size_t x_steps,
    std::size_t t_steps,
    double wave_speed = 1.0);

Grid2D poisson_dirichlet(
    const std::function<double(double, double)>& source,
    const std::function<double(double, double)>& boundary,
    double x_left,
    double x_right,
    double y_bottom,
    double y_top,
    std::size_t x_interior,
    std::size_t y_interior);

} // namespace numerical
