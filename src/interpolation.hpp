#pragma once

#include "types.hpp"

namespace numerical {

double lagrange_interpolate(
    const Vector& xs,
    const Vector& ys,
    double x);

double neville_interpolate(
    const Vector& xs,
    const Vector& ys,
    double x);

Vector divided_difference_coefficients(
    const Vector& xs,
    const Vector& ys);

double evaluate_newton_polynomial(
    const Vector& xs,
    const Vector& coefficients,
    double x);

double hermite_interpolate(
    const Vector& xs,
    const Vector& ys,
    const Vector& derivatives,
    double x);

CubicSpline natural_cubic_spline(
    const Vector& xs,
    const Vector& ys);

CubicSpline clamped_cubic_spline(
    const Vector& xs,
    const Vector& ys,
    double left_derivative,
    double right_derivative);

double evaluate_spline(const CubicSpline& spline, double x);

} // namespace numerical
