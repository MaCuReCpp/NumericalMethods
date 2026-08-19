#include "differentiation_integration.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>

namespace numerical {

namespace {

void require_nonzero_step(double step)
{
    if (step <= 0.0) {
        throw std::invalid_argument("step must be positive");
    }
}

} // namespace

double finite_difference(
    const ScalarFunction& f,
    double x,
    double step,
    DifferenceMethod method)
{
    require_nonzero_step(step);

    switch (method) {
    case DifferenceMethod::Forward:
        return (f(x + step) - f(x)) / step;
    case DifferenceMethod::Backward:
        return (f(x) - f(x - step)) / step;
    case DifferenceMethod::Central:
        return (f(x + step) - f(x - step)) / (2.0 * step);
    }

    throw std::invalid_argument("unknown difference method");
}

double trapezoidal_rule(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t subintervals)
{
    if (subintervals == 0) {
        throw std::invalid_argument("subintervals must be positive");
    }

    const double h = (b - a) / static_cast<double>(subintervals);
    double sum = 0.5 * (f(a) + f(b));

    for (std::size_t i = 1; i < subintervals; ++i) {
        sum += f(a + h * static_cast<double>(i));
    }

    return h * sum;
}

double simpson_rule(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t subintervals)
{
    if (subintervals == 0 || subintervals % 2 != 0) {
        throw std::invalid_argument("Simpson rule requires a positive even number of subintervals");
    }

    const double h = (b - a) / static_cast<double>(subintervals);
    double sum = f(a) + f(b);

    for (std::size_t i = 1; i < subintervals; ++i) {
        const double weight = (i % 2 == 0) ? 2.0 : 4.0;
        sum += weight * f(a + h * static_cast<double>(i));
    }

    return (h / 3.0) * sum;
}


double adaptive_simpson_recursive(
    const ScalarFunction& f,
    double a,
    double b,
    double fa,
    double fm,
    double fb,
    double whole,
    double tolerance,
    std::size_t depth)
{
    const double mid = 0.5 * (a + b);
    const double left_mid = 0.5 * (a + mid);
    const double right_mid = 0.5 * (mid + b);
    const double f_left_mid = f(left_mid);
    const double f_right_mid = f(right_mid);

    const double left = (mid - a) * (fa + 4.0 * f_left_mid + fm) / 6.0;
    const double right = (b - mid) * (fm + 4.0 * f_right_mid + fb) / 6.0;
    const double refined = left + right;

    if (depth == 0 || std::abs(refined - whole) <= 15.0 * tolerance) {
        return refined + (refined - whole) / 15.0;
    }

    return adaptive_simpson_recursive(f, a, mid, fa, f_left_mid, fm, left, 0.5 * tolerance, depth - 1) +
           adaptive_simpson_recursive(f, mid, b, fm, f_right_mid, fb, right, 0.5 * tolerance, depth - 1);
}


double richardson_derivative(
    const ScalarFunction& f,
    double x,
    double initial_step,
    std::size_t levels)
{
    if (initial_step <= 0.0 || levels == 0) {
        throw std::invalid_argument("step and levels must be positive");
    }

    Matrix table(levels, Vector(levels, 0.0));
    for (std::size_t i = 0; i < levels; ++i) {
        const double h = initial_step / static_cast<double>(1ULL << i);
        table[i][0] = finite_difference(f, x, h, DifferenceMethod::Central);
    }

    for (std::size_t j = 1; j < levels; ++j) {
        const double factor = std::pow(4.0, static_cast<double>(j));
        for (std::size_t i = j; i < levels; ++i) {
            table[i][j] = table[i][j - 1] + (table[i][j - 1] - table[i - 1][j - 1]) / (factor - 1.0);
        }
    }

    return table.back().back();
}

double romberg_integration(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t levels)
{
    if (levels == 0) {
        throw std::invalid_argument("levels must be positive");
    }

    Matrix table(levels, Vector(levels, 0.0));
    for (std::size_t i = 0; i < levels; ++i) {
        table[i][0] = trapezoidal_rule(f, a, b, static_cast<std::size_t>(1ULL << i));
    }

    for (std::size_t j = 1; j < levels; ++j) {
        const double factor = std::pow(4.0, static_cast<double>(j));
        for (std::size_t i = j; i < levels; ++i) {
            table[i][j] = table[i][j - 1] + (table[i][j - 1] - table[i - 1][j - 1]) / (factor - 1.0);
        }
    }

    return table.back().back();
}

double adaptive_simpson(
    const ScalarFunction& f,
    double a,
    double b,
    double tolerance,
    std::size_t max_depth)
{
    if (tolerance <= 0.0) {
        throw std::invalid_argument("tolerance must be positive");
    }
    const double mid = 0.5 * (a + b);
    const double fa = f(a);
    const double fm = f(mid);
    const double fb = f(b);
    const double whole = (b - a) * (fa + 4.0 * fm + fb) / 6.0;
    return adaptive_simpson_recursive(f, a, b, fa, fm, fb, whole, tolerance, max_depth);
}

double gaussian_quadrature(
    const ScalarFunction& f,
    double a,
    double b,
    std::size_t points)
{
    static const Vector nodes2 = {-0.5773502691896257, 0.5773502691896257};
    static const Vector weights2 = {1.0, 1.0};
    static const Vector nodes3 = {-0.7745966692414834, 0.0, 0.7745966692414834};
    static const Vector weights3 = {0.5555555555555556, 0.8888888888888888, 0.5555555555555556};
    static const Vector nodes4 = {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526};
    static const Vector weights4 = {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538};
    static const Vector nodes5 = {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640};
    static const Vector weights5 = {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891};

    const Vector* nodes = nullptr;
    const Vector* weights = nullptr;
    if (points == 2) {
        nodes = &nodes2;
        weights = &weights2;
    } else if (points == 3) {
        nodes = &nodes3;
        weights = &weights3;
    } else if (points == 4) {
        nodes = &nodes4;
        weights = &weights4;
    } else if (points == 5) {
        nodes = &nodes5;
        weights = &weights5;
    } else {
        throw std::invalid_argument("Gaussian quadrature supports 2 to 5 points");
    }

    const double center = 0.5 * (a + b);
    const double half_width = 0.5 * (b - a);
    double sum = 0.0;
    for (std::size_t i = 0; i < points; ++i) {
        sum += (*weights)[i] * f(center + half_width * (*nodes)[i]);
    }
    return half_width * sum;
}

} // namespace numerical
