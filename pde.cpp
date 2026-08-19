#include "pde.hpp"
#include "linear_algebra.hpp"

#include <cmath>

namespace numerical {

Grid1D finite_difference_bvp(
    const ScalarFunction& p,
    const ScalarFunction& q,
    const ScalarFunction& r,
    double a,
    double b,
    double alpha,
    double beta,
    std::size_t interior_points)
{
    if (interior_points == 0) {
        throw std::invalid_argument("interior_points must be positive");
    }
    const std::size_t n = interior_points;
    const double h = (b - a) / static_cast<double>(n + 1);
    Matrix mat(n, Vector(n, 0.0));
    Vector rhs(n, 0.0);

    for (std::size_t i = 0; i < n; ++i) {
        const double x = a + h * static_cast<double>(i + 1);
        const double lower = 1.0 + 0.5 * h * p(x);
        const double diag = -2.0 - h * h * q(x);
        const double upper = 1.0 - 0.5 * h * p(x);
        if (i > 0) {
            mat[i][i - 1] = lower;
        } else {
            rhs[i] -= lower * alpha;
        }
        mat[i][i] = diag;
        if (i + 1 < n) {
            mat[i][i + 1] = upper;
        } else {
            rhs[i] -= upper * beta;
        }
        rhs[i] += h * h * r(x);
    }

    const Vector interior = gaussian_elimination(mat, rhs).solution;
    Grid1D grid;
    grid.points.resize(n + 2);
    grid.values.resize(n + 2);
    grid.points.front() = a;
    grid.values.front() = alpha;
    for (std::size_t i = 0; i < n; ++i) {
        grid.points[i + 1] = a + h * static_cast<double>(i + 1);
        grid.values[i + 1] = interior[i];
    }
    grid.points.back() = b;
    grid.values.back() = beta;
    return grid;
}

HeatEquationResult heat_equation_explicit(
    const ScalarFunction& initial,
    const ScalarFunction& left_boundary,
    const ScalarFunction& right_boundary,
    double x_left,
    double x_right,
    double t_final,
    std::size_t x_steps,
    std::size_t t_steps,
    double diffusivity)
{
    if (x_steps < 2 || t_steps == 0 || diffusivity <= 0.0) {
        throw std::invalid_argument("invalid heat equation grid or diffusivity");
    }
    const double dx = (x_right - x_left) / static_cast<double>(x_steps);
    const double dt = t_final / static_cast<double>(t_steps);
    const double lambda = diffusivity * dt / (dx * dx);
    if (lambda > 0.5) {
        throw std::invalid_argument("explicit heat equation scheme requires diffusivity * dt / dx^2 <= 0.5");
    }

    HeatEquationResult result;
    result.x.resize(x_steps + 1);
    result.t.resize(t_steps + 1);
    result.values.assign(t_steps + 1, Vector(x_steps + 1, 0.0));

    for (std::size_t i = 0; i <= x_steps; ++i) {
        result.x[i] = x_left + dx * static_cast<double>(i);
        result.values[0][i] = initial(result.x[i]);
    }
    for (std::size_t n = 0; n <= t_steps; ++n) {
        result.t[n] = dt * static_cast<double>(n);
        result.values[n][0] = left_boundary(result.t[n]);
        result.values[n][x_steps] = right_boundary(result.t[n]);
    }
    for (std::size_t n = 0; n < t_steps; ++n) {
        for (std::size_t i = 1; i < x_steps; ++i) {
            result.values[n + 1][i] = result.values[n][i] +
                lambda * (result.values[n][i - 1] - 2.0 * result.values[n][i] + result.values[n][i + 1]);
        }
    }
    return result;
}

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
    double wave_speed)
{
    if (x_steps < 2 || t_steps == 0 || wave_speed <= 0.0) {
        throw std::invalid_argument("invalid wave equation grid or speed");
    }
    const double dx = (x_right - x_left) / static_cast<double>(x_steps);
    const double dt = t_final / static_cast<double>(t_steps);
    const double lambda = wave_speed * dt / dx;
    if (lambda > 1.0) {
        throw std::invalid_argument("explicit wave scheme requires wave_speed * dt / dx <= 1");
    }

    WaveEquationResult result;
    result.x.resize(x_steps + 1);
    result.t.resize(t_steps + 1);
    result.values.assign(t_steps + 1, Vector(x_steps + 1, 0.0));
    const double lambda2 = lambda * lambda;

    for (std::size_t i = 0; i <= x_steps; ++i) {
        result.x[i] = x_left + dx * static_cast<double>(i);
        result.values[0][i] = initial_displacement(result.x[i]);
    }
    for (std::size_t n = 0; n <= t_steps; ++n) {
        result.t[n] = dt * static_cast<double>(n);
        result.values[n][0] = left_boundary(result.t[n]);
        result.values[n][x_steps] = right_boundary(result.t[n]);
    }
    for (std::size_t i = 1; i < x_steps; ++i) {
        result.values[1][i] = result.values[0][i] + dt * initial_velocity(result.x[i]) +
            0.5 * lambda2 * (result.values[0][i - 1] - 2.0 * result.values[0][i] + result.values[0][i + 1]);
    }
    for (std::size_t n = 1; n < t_steps; ++n) {
        for (std::size_t i = 1; i < x_steps; ++i) {
            result.values[n + 1][i] = 2.0 * (1.0 - lambda2) * result.values[n][i] -
                result.values[n - 1][i] + lambda2 * (result.values[n][i - 1] + result.values[n][i + 1]);
        }
    }
    return result;
}

Grid2D poisson_dirichlet(
    const std::function<double(double, double)>& source,
    const std::function<double(double, double)>& boundary,
    double x_left,
    double x_right,
    double y_bottom,
    double y_top,
    std::size_t x_interior,
    std::size_t y_interior)
{
    if (x_interior == 0 || y_interior == 0) {
        throw std::invalid_argument("interior grid sizes must be positive");
    }
    const double hx = (x_right - x_left) / static_cast<double>(x_interior + 1);
    const double hy = (y_top - y_bottom) / static_cast<double>(y_interior + 1);
    if (std::abs(hx - hy) > 1e-12) {
        throw std::invalid_argument("poisson_dirichlet currently requires equal grid spacing in x and y");
    }
    const std::size_t unknowns = x_interior * y_interior;
    Matrix mat(unknowns, Vector(unknowns, 0.0));
    Vector rhs(unknowns, 0.0);
    auto index = [=](std::size_t i, std::size_t j) {
        return j * x_interior + i;
    };

    for (std::size_t j = 0; j < y_interior; ++j) {
        for (std::size_t i = 0; i < x_interior; ++i) {
            const std::size_t row = index(i, j);
            const double x = x_left + hx * static_cast<double>(i + 1);
            const double y = y_bottom + hy * static_cast<double>(j + 1);
            mat[row][row] = -4.0;
            rhs[row] = hx * hx * source(x, y);

            const double neighbor_x[4] = {x - hx, x + hx, x, x};
            const double neighbor_y[4] = {y, y, y - hy, y + hy};
            const int di[4] = {-1, 1, 0, 0};
            const int dj[4] = {0, 0, -1, 1};
            for (std::size_t side = 0; side < 4; ++side) {
                const int ni = static_cast<int>(i) + di[side];
                const int nj = static_cast<int>(j) + dj[side];
                if (ni >= 0 && ni < static_cast<int>(x_interior) && nj >= 0 && nj < static_cast<int>(y_interior)) {
                    mat[row][index(static_cast<std::size_t>(ni), static_cast<std::size_t>(nj))] = 1.0;
                } else {
                    rhs[row] -= boundary(neighbor_x[side], neighbor_y[side]);
                }
            }
        }
    }

    const Vector interior = gaussian_elimination(mat, rhs).solution;
    Grid2D grid;
    grid.x.resize(x_interior + 2);
    grid.y.resize(y_interior + 2);
    grid.values.assign(y_interior + 2, Vector(x_interior + 2, 0.0));
    for (std::size_t i = 0; i < grid.x.size(); ++i) {
        grid.x[i] = x_left + hx * static_cast<double>(i);
    }
    for (std::size_t j = 0; j < grid.y.size(); ++j) {
        grid.y[j] = y_bottom + hy * static_cast<double>(j);
    }
    for (std::size_t j = 0; j < grid.y.size(); ++j) {
        for (std::size_t i = 0; i < grid.x.size(); ++i) {
            const bool is_boundary = i == 0 || j == 0 || i + 1 == grid.x.size() || j + 1 == grid.y.size();
            grid.values[j][i] = is_boundary ? boundary(grid.x[i], grid.y[j]) : interior[index(i - 1, j - 1)];
        }
    }
    return grid;
}

} // namespace numerical
