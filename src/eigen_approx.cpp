#include "eigen_approx.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>

namespace numerical {
using detail::dot;
using detail::multiply;
using detail::require_square;
using detail::subtract;

Matrix transpose_matrix(const Matrix& a)
{
    if (a.empty()) {
        return {};
    }
    Matrix result(a.front().size(), Vector(a.size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < a[i].size(); ++j) {
            result[j][i] = a[i][j];
        }
    }
    return result;
}

Matrix multiply_matrix(const Matrix& a, const Matrix& b)
{
    if (a.empty() || b.empty() || a.front().size() != b.size()) {
        throw std::invalid_argument("matrix dimensions are incompatible");
    }
    Matrix result(a.size(), Vector(b.front().size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t k = 0; k < b.size(); ++k) {
            for (std::size_t j = 0; j < b.front().size(); ++j) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

Matrix identity_matrix(std::size_t n)
{
    Matrix result(n, Vector(n, 0.0));
    for (std::size_t i = 0; i < n; ++i) {
        result[i][i] = 1.0;
    }
    return result;
}

double off_diagonal_norm(const Matrix& a)
{
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < a.size(); ++j) {
            if (i != j) {
                sum += a[i][j] * a[i][j];
            }
        }
    }
    return std::sqrt(sum);
}


EigenResult power_method(
    const Matrix& a,
    Vector initial_guess,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    require_square(a);
    if (initial_guess.size() != a.size()) {
        throw std::invalid_argument("initial vector dimension is incompatible");
    }

    double norm = l2_norm(initial_guess);
    if (norm < detail::kTiny) {
        throw std::invalid_argument("initial vector must be non-zero");
    }
    for (double& value : initial_guess) {
        value /= norm;
    }

    Vector x = initial_guess;
    double eigenvalue = 0.0;
    for (std::size_t iteration = 1; iteration <= max_iterations; ++iteration) {
        Vector y = multiply(a, x);
        norm = l2_norm(y);
        if (norm < detail::kTiny) {
            throw std::runtime_error("power method produced a near-zero vector");
        }
        for (double& value : y) {
            value /= norm;
        }
        const Vector ay = multiply(a, y);
        const double next_eigenvalue = dot(y, ay);
        Vector residual = ay;
        for (std::size_t i = 0; i < residual.size(); ++i) {
            residual[i] -= next_eigenvalue * y[i];
        }
        const double residual_norm = l2_norm(residual);
        if (std::abs(next_eigenvalue - eigenvalue) <= tolerance.absolute + tolerance.relative * std::max(1.0, std::abs(next_eigenvalue)) ||
            residual_norm <= tolerance.absolute) {
            return {next_eigenvalue, y, residual_norm, iteration, true};
        }
        eigenvalue = next_eigenvalue;
        x = y;
    }
    return {eigenvalue, x, l2_norm(subtract(multiply(a, x), Vector(a.size(), 0.0))), max_iterations, false};
}

QRDecomposition qr_decomposition(const Matrix& a)
{
    if (a.empty()) {
        throw std::invalid_argument("matrix must be non-empty");
    }
    const std::size_t rows = a.size();
    const std::size_t cols = a.front().size();
    Matrix q(rows, Vector(cols, 0.0));
    Matrix r(cols, Vector(cols, 0.0));

    for (std::size_t k = 0; k < cols; ++k) {
        Vector v(rows, 0.0);
        for (std::size_t i = 0; i < rows; ++i) {
            v[i] = a[i][k];
        }
        for (std::size_t j = 0; j < k; ++j) {
            Vector qj(rows, 0.0);
            for (std::size_t i = 0; i < rows; ++i) {
                qj[i] = q[i][j];
            }
            r[j][k] = dot(qj, v);
            for (std::size_t i = 0; i < rows; ++i) {
                v[i] -= r[j][k] * q[i][j];
            }
        }
        r[k][k] = l2_norm(v);
        if (r[k][k] < detail::kTiny) {
            throw std::runtime_error("QR decomposition encountered dependent columns");
        }
        for (std::size_t i = 0; i < rows; ++i) {
            q[i][k] = v[i] / r[k][k];
        }
    }

    return {q, r};
}

Vector qr_eigenvalues(
    Matrix a,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    require_square(a);
    for (std::size_t iteration = 0; iteration < max_iterations; ++iteration) {
        const auto qr = qr_decomposition(a);
        a = multiply_matrix(qr.r, qr.q);
        if (off_diagonal_norm(a) <= tolerance.absolute) {
            break;
        }
    }
    Vector values(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        values[i] = a[i][i];
    }
    return values;
}

SVDResult svd(
    const Matrix& a,
    Tolerance tolerance,
    std::size_t max_iterations)
{
    if (a.empty()) {
        throw std::invalid_argument("matrix must be non-empty");
    }
    const std::size_t rows = a.size();
    const std::size_t cols = a.front().size();
    Matrix ata = multiply_matrix(transpose_matrix(a), a);
    Matrix v(cols, Vector(cols, 0.0));
    Vector singular_values(cols, 0.0);

    for (std::size_t k = 0; k < cols; ++k) {
        Vector guess(cols, 1.0);
        guess[k] = 2.0;
        const auto eigen = power_method(ata, guess, tolerance, max_iterations);
        const double lambda = std::max(0.0, eigen.eigenvalue);
        singular_values[k] = std::sqrt(lambda);
        for (std::size_t i = 0; i < cols; ++i) {
            v[i][k] = eigen.eigenvector[i];
        }
        for (std::size_t i = 0; i < cols; ++i) {
            for (std::size_t j = 0; j < cols; ++j) {
                ata[i][j] -= lambda * eigen.eigenvector[i] * eigen.eigenvector[j];
            }
        }
    }

    Matrix u(rows, Vector(cols, 0.0));
    for (std::size_t k = 0; k < cols; ++k) {
        if (singular_values[k] <= detail::kTiny) {
            continue;
        }
        Vector vk(cols, 0.0);
        for (std::size_t i = 0; i < cols; ++i) {
            vk[i] = v[i][k];
        }
        Vector av = multiply(a, vk);
        for (std::size_t i = 0; i < rows; ++i) {
            u[i][k] = av[i] / singular_values[k];
        }
    }

    return {u, singular_values, transpose_matrix(v)};
}

} // namespace numerical
