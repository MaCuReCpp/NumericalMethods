#include "differentiation_integration.hpp"
#include "errors.hpp"
#include "interpolation.hpp"
#include "linear_algebra.hpp"
#include "ode.hpp"
#include "roots.hpp"
#include "special_functions.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace numerical {
namespace {

constexpr double kPivotEpsilon = 1e-14;
constexpr double kSeriesTolerance = 1e-16;
constexpr std::size_t kMaxSeriesTerms = 200;

void require_nonzero_step(double step)
{
    if (step <= 0.0) {
        throw std::invalid_argument("step must be positive");
    }
}

void require_same_size(const Vector& xs, const Vector& ys)
{
    if (xs.size() != ys.size() || xs.empty()) {
        throw std::invalid_argument("vectors must be non-empty and have the same size");
    }
}

void require_square_matrix(const Matrix& a, const Vector& b)
{
    if (a.empty() || a.size() != b.size()) {
        throw std::invalid_argument("matrix and vector dimensions are incompatible");
    }

    for (const auto& row : a) {
        if (row.size() != a.size()) {
            throw std::invalid_argument("matrix must be square");
        }
    }
}

bool close_enough(double current, double previous, double residual, Tolerance tolerance)
{
    const double scale = std::max(1.0, std::abs(current));
    const double step_error = std::abs(current - previous);
    return std::abs(residual) <= tolerance.absolute ||
           step_error <= tolerance.absolute + tolerance.relative * scale;
}

Vector matvec(const Matrix& a, const Vector& x)
{
    Vector result(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < x.size(); ++j) {
            result[i] += a[i][j] * x[j];
        }
    }
    return result;
}

Vector residual_vector(const Matrix& a, const Vector& x, const Vector& b)
{
    Vector ax = matvec(a, x);
    for (std::size_t i = 0; i < ax.size(); ++i) {
        ax[i] -= b[i];
    }
    return ax;
}

double bessel_j0_series(double x)
{
    const double half_x_squared = 0.25 * x * x;
    double term = 1.0;
    double sum = term;

    for (std::size_t k = 1; k <= kMaxSeriesTerms; ++k) {
        const double denominator = static_cast<double>(k * k);
        term *= -half_x_squared / denominator;
        sum += term;

        if (std::abs(term) <= kSeriesTolerance * std::max(1.0, std::abs(sum))) {
            break;
        }
    }

    return sum;
}

double bessel_j1_series(double x)
{
    const double half_x_squared = 0.25 * x * x;
    double term = 0.5 * x;
    double sum = term;

    for (std::size_t k = 1; k <= kMaxSeriesTerms; ++k) {
        const double denominator = static_cast<double>(k * (k + 1));
        term *= -half_x_squared / denominator;
        sum += term;

        if (std::abs(term) <= kSeriesTolerance * std::max(1.0, std::abs(sum))) {
            break;
        }
    }

    return sum;
}

} // namespace

double absolute_error(double exact, double approximate)
{
    return std::abs(exact - approximate);
}

double relative_error(double exact, double approximate)
{
    if (exact == 0.0) {
        throw std::invalid_argument("relative error is undefined when exact value is zero");
    }
    return std::abs((exact - approximate) / exact);
}

double l2_norm(const Vector& values)
{
    double sum = 0.0;
    for (double value : values) {
        sum += value * value;
    }
    return std::sqrt(sum);
}

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

double bessel_j0(double x)
{
    return bessel_j0_series(x);
}

double bessel_j1(double x)
{
    return bessel_j1_series(x);
}

double bessel_j(int order, double x)
{
    if (order < 0) {
        const int positive_order = -order;
        const double value = bessel_j(positive_order, x);
        return (positive_order % 2 == 0) ? value : -value;
    }

    if (order == 0) {
        return bessel_j0(x);
    }
    if (order == 1) {
        return bessel_j1(x);
    }
    if (x == 0.0) {
        return 0.0;
    }

    double previous = bessel_j0(x);
    double current = bessel_j1(x);
    for (int n = 1; n < order; ++n) {
        const double next = (2.0 * static_cast<double>(n) / x) * current - previous;
        previous = current;
        current = next;
    }

    return current;
}

double lagrange_interpolate(
    const Vector& xs,
    const Vector& ys,
    double x)
{
    require_same_size(xs, ys);

    double result = 0.0;
    for (std::size_t i = 0; i < xs.size(); ++i) {
        double basis = 1.0;
        for (std::size_t j = 0; j < xs.size(); ++j) {
            if (i == j) {
                continue;
            }
            const double denominator = xs[i] - xs[j];
            if (std::abs(denominator) < kPivotEpsilon) {
                throw std::invalid_argument("interpolation nodes must be distinct");
            }
            basis *= (x - xs[j]) / denominator;
        }
        result += ys[i] * basis;
    }

    return result;
}

Vector divided_difference_coefficients(
    const Vector& xs,
    const Vector& ys)
{
    require_same_size(xs, ys);

    Vector coefficients = ys;
    for (std::size_t order = 1; order < xs.size(); ++order) {
        for (std::size_t i = xs.size() - 1; i >= order; --i) {
            const double denominator = xs[i] - xs[i - order];
            if (std::abs(denominator) < kPivotEpsilon) {
                throw std::invalid_argument("interpolation nodes must be distinct");
            }
            coefficients[i] = (coefficients[i] - coefficients[i - 1]) / denominator;
        }
    }

    return coefficients;
}

double evaluate_newton_polynomial(
    const Vector& xs,
    const Vector& coefficients,
    double x)
{
    if (xs.empty() || xs.size() != coefficients.size()) {
        throw std::invalid_argument("nodes and coefficients must have the same non-empty size");
    }

    double result = coefficients.back();
    for (std::size_t i = coefficients.size() - 1; i-- > 0;) {
        result = result * (x - xs[i]) + coefficients[i];
    }

    return result;
}

LinearSolveResult gaussian_elimination(Matrix a, Vector b)
{
    require_square_matrix(a, b);

    const Matrix original_a = a;
    const Vector original_b = b;
    const std::size_t n = a.size();

    for (std::size_t pivot = 0; pivot < n; ++pivot) {
        std::size_t max_row = pivot;
        for (std::size_t row = pivot + 1; row < n; ++row) {
            if (std::abs(a[row][pivot]) > std::abs(a[max_row][pivot])) {
                max_row = row;
            }
        }

        if (std::abs(a[max_row][pivot]) < kPivotEpsilon) {
            throw std::runtime_error("matrix is singular or ill-conditioned");
        }

        std::swap(a[pivot], a[max_row]);
        std::swap(b[pivot], b[max_row]);

        for (std::size_t row = pivot + 1; row < n; ++row) {
            const double factor = a[row][pivot] / a[pivot][pivot];
            a[row][pivot] = 0.0;
            for (std::size_t col = pivot + 1; col < n; ++col) {
                a[row][col] -= factor * a[pivot][col];
            }
            b[row] -= factor * b[pivot];
        }
    }

    Vector x(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double sum = b[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= a[i][j] * x[j];
        }
        x[i] = sum / a[i][i];
    }

    return {x, l2_norm(residual_vector(original_a, x, original_b)), n, true};
}

LinearSolveResult jacobi(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    require_square_matrix(a, b);
    if (initial_guess.size() != b.size()) {
        throw std::invalid_argument("initial guess dimension is incompatible");
    }

    const std::size_t n = a.size();
    Vector x = initial_guess;
    Vector next = x;

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        for (std::size_t i = 0; i < n; ++i) {
            if (std::abs(a[i][i]) < kPivotEpsilon) {
                throw std::runtime_error("Jacobi method encountered a near-zero diagonal entry");
            }

            double sum = b[i];
            for (std::size_t j = 0; j < n; ++j) {
                if (i != j) {
                    sum -= a[i][j] * x[j];
                }
            }
            next[i] = sum / a[i][i];
        }

        Vector difference(n, 0.0);
        for (std::size_t i = 0; i < n; ++i) {
            difference[i] = next[i] - x[i];
        }

        const double residual = l2_norm(residual_vector(a, next, b));
        if (residual <= tolerance.absolute ||
            l2_norm(difference) <= tolerance.absolute + tolerance.relative * std::max(1.0, l2_norm(next))) {
            return {next, residual, iteration, true};
        }

        x = next;
    }

    return {x, l2_norm(residual_vector(a, x, b)), max_iterations, false};
}

std::vector<OdePoint> euler_method(
    const OdeFunction& f,
    double t0,
    double y0,
    double step,
    std::size_t steps)
{
    require_nonzero_step(step);

    std::vector<OdePoint> solution;
    solution.reserve(steps + 1);
    solution.push_back({t0, y0});

    double t = t0;
    double y = y0;
    for (std::size_t i = 0; i < steps; ++i) {
        y += step * f(t, y);
        t += step;
        solution.push_back({t, y});
    }

    return solution;
}

std::vector<OdePoint> runge_kutta4(
    const OdeFunction& f,
    double t0,
    double y0,
    double step,
    std::size_t steps)
{
    require_nonzero_step(step);

    std::vector<OdePoint> solution;
    solution.reserve(steps + 1);
    solution.push_back({t0, y0});

    double t = t0;
    double y = y0;
    for (std::size_t i = 0; i < steps; ++i) {
        const double k1 = f(t, y);
        const double k2 = f(t + 0.5 * step, y + 0.5 * step * k1);
        const double k3 = f(t + 0.5 * step, y + 0.5 * step * k2);
        const double k4 = f(t + step, y + step * k3);

        y += (step / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
        t += step;
        solution.push_back({t, y});
    }

    return solution;
}

} // namespace numerical
