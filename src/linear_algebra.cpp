#include "linear_algebra.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>

namespace numerical {
using detail::dot;
using detail::l2_distance;
using detail::multiply;
using detail::require_square;
using detail::subtract;

LUDecomposition lu_decomposition(Matrix a)
{
    require_square(a);
    const std::size_t n = a.size();
    Matrix lower(n, Vector(n, 0.0));
    std::vector<std::size_t> permutation(n, 0);
    for (std::size_t i = 0; i < n; ++i) {
        lower[i][i] = 1.0;
        permutation[i] = i;
    }

    int swaps = 0;
    for (std::size_t k = 0; k < n; ++k) {
        std::size_t pivot = k;
        for (std::size_t i = k + 1; i < n; ++i) {
            if (std::abs(a[i][k]) > std::abs(a[pivot][k])) {
                pivot = i;
            }
        }
        if (std::abs(a[pivot][k]) < detail::kTiny) {
            throw std::runtime_error("matrix is singular or ill-conditioned");
        }
        if (pivot != k) {
            std::swap(a[pivot], a[k]);
            std::swap(permutation[pivot], permutation[k]);
            for (std::size_t j = 0; j < k; ++j) {
                std::swap(lower[pivot][j], lower[k][j]);
            }
            ++swaps;
        }

        for (std::size_t i = k + 1; i < n; ++i) {
            lower[i][k] = a[i][k] / a[k][k];
            for (std::size_t j = k; j < n; ++j) {
                a[i][j] -= lower[i][k] * a[k][j];
            }
        }
    }

    Matrix upper(n, Vector(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = i; j < n; ++j) {
            upper[i][j] = a[i][j];
        }
    }

    return {lower, upper, permutation, swaps};
}

Vector solve_lu(const LUDecomposition& lu, const Vector& b)
{
    const std::size_t n = lu.lower.size();
    if (b.size() != n || lu.upper.size() != n || lu.permutation.size() != n) {
        throw std::invalid_argument("LU decomposition and vector dimensions are incompatible");
    }

    Vector pb(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        pb[i] = b[lu.permutation[i]];
    }

    Vector y(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        double sum = pb[i];
        for (std::size_t j = 0; j < i; ++j) {
            sum -= lu.lower[i][j] * y[j];
        }
        y[i] = sum;
    }

    Vector x(n, 0.0);
    for (std::size_t i = n; i-- > 0;) {
        double sum = y[i];
        for (std::size_t j = i + 1; j < n; ++j) {
            sum -= lu.upper[i][j] * x[j];
        }
        if (std::abs(lu.upper[i][i]) < detail::kTiny) {
            throw std::runtime_error("LU decomposition has a near-zero pivot");
        }
        x[i] = sum / lu.upper[i][i];
    }
    return x;
}

double determinant(Matrix a)
{
    const auto lu = lu_decomposition(a);
    double det = (lu.row_swaps % 2 == 0) ? 1.0 : -1.0;
    for (std::size_t i = 0; i < lu.upper.size(); ++i) {
        det *= lu.upper[i][i];
    }
    return det;
}

Matrix inverse(Matrix a)
{
    require_square(a);
    const std::size_t n = a.size();
    const auto lu = lu_decomposition(a);
    Matrix inv(n, Vector(n, 0.0));
    for (std::size_t col = 0; col < n; ++col) {
        Vector e(n, 0.0);
        e[col] = 1.0;
        const Vector solution = solve_lu(lu, e);
        for (std::size_t row = 0; row < n; ++row) {
            inv[row][col] = solution[row];
        }
    }
    return inv;
}

LinearSolveResult gauss_seidel(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    require_square(a);
    if (b.size() != a.size() || initial_guess.size() != a.size()) {
        throw std::invalid_argument("dimensions are incompatible");
    }

    Vector x = initial_guess;
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const Vector previous = x;
        for (std::size_t i = 0; i < a.size(); ++i) {
            double sum = b[i];
            for (std::size_t j = 0; j < a.size(); ++j) {
                if (i != j) {
                    sum -= a[i][j] * x[j];
                }
            }
            if (std::abs(a[i][i]) < detail::kTiny) {
                throw std::runtime_error("Gauss-Seidel method encountered a near-zero diagonal entry");
            }
            x[i] = sum / a[i][i];
        }
        const double residual = l2_norm(subtract(multiply(a, x), b));
        if (residual <= tolerance.absolute ||
            l2_distance(x, previous) <= tolerance.absolute + tolerance.relative * std::max(1.0, l2_norm(x))) {
            return {x, residual, iteration, true};
        }
    }
    return {x, l2_norm(subtract(multiply(a, x), b)), max_iterations, false};
}

LinearSolveResult sor(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    double omega,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    if (omega <= 0.0 || omega >= 2.0) {
        throw std::invalid_argument("omega should be in (0, 2)");
    }
    require_square(a);
    if (b.size() != a.size() || initial_guess.size() != a.size()) {
        throw std::invalid_argument("dimensions are incompatible");
    }

    Vector x = initial_guess;
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const Vector previous = x;
        for (std::size_t i = 0; i < a.size(); ++i) {
            double sum = b[i];
            for (std::size_t j = 0; j < a.size(); ++j) {
                if (i != j) {
                    sum -= a[i][j] * x[j];
                }
            }
            if (std::abs(a[i][i]) < detail::kTiny) {
                throw std::runtime_error("SOR method encountered a near-zero diagonal entry");
            }
            x[i] = (1.0 - omega) * x[i] + omega * sum / a[i][i];
        }
        const double residual = l2_norm(subtract(multiply(a, x), b));
        if (residual <= tolerance.absolute ||
            l2_distance(x, previous) <= tolerance.absolute + tolerance.relative * std::max(1.0, l2_norm(x))) {
            return {x, residual, iteration, true};
        }
    }
    return {x, l2_norm(subtract(multiply(a, x), b)), max_iterations, false};
}


LinearSolveResult conjugate_gradient(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    require_square(a);
    if (b.size() != a.size() || initial_guess.size() != a.size()) {
        throw std::invalid_argument("dimensions are incompatible");
    }

    Vector x = initial_guess;
    Vector r = subtract(b, multiply(a, x));
    Vector p = r;
    double rs_old = dot(r, r);

    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const Vector ap = multiply(a, p);
        const double denominator = dot(p, ap);
        if (std::abs(denominator) < detail::kTiny) {
            throw std::runtime_error("conjugate gradient encountered a near-zero denominator");
        }
        const double alpha = rs_old / denominator;
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] += alpha * p[i];
            r[i] -= alpha * ap[i];
        }
        const double rs_new = dot(r, r);
        const double residual = std::sqrt(rs_new);
        if (residual <= tolerance.absolute) {
            return {x, residual, iteration, true};
        }
        const double beta = rs_new / rs_old;
        for (std::size_t i = 0; i < p.size(); ++i) {
            p[i] = r[i] + beta * p[i];
        }
        rs_old = rs_new;
    }

    return {x, std::sqrt(rs_old), max_iterations, false};
}

LinearSolveResult iterative_refinement(
    const Matrix& a,
    const Vector& b,
    Vector initial_solution,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    require_square(a);
    if (b.size() != a.size() || initial_solution.size() != a.size()) {
        throw std::invalid_argument("dimensions are incompatible");
    }

    const auto lu = lu_decomposition(a);
    Vector x = initial_solution;
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        const Vector residual_vec = subtract(b, multiply(a, x));
        const double residual = l2_norm(residual_vec);
        if (residual <= tolerance.absolute) {
            return {x, residual, iteration - 1, true};
        }
        const Vector correction = solve_lu(lu, residual_vec);
        for (std::size_t i = 0; i < x.size(); ++i) {
            x[i] += correction[i];
        }
        if (l2_norm(correction) <= tolerance.absolute + tolerance.relative * std::max(1.0, l2_norm(x))) {
            return {x, l2_norm(subtract(b, multiply(a, x))), iteration, true};
        }
    }
    return {x, l2_norm(subtract(b, multiply(a, x))), max_iterations, false};
}

} // namespace numerical
