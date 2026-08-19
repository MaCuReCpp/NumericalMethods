#include "approximation.hpp"
#include "internal.hpp"
#include "linear_algebra.hpp"

#include <algorithm>
#include <cmath>

namespace numerical {
using detail::require_same_size_nonempty;

namespace {
constexpr double kPi = 3.141592653589793238462643383279502884;
}

Vector least_squares_polynomial(
    const Vector& xs,
    const Vector& ys,
    std::size_t degree)
{
    require_same_size_nonempty(xs, ys);
    const std::size_t m = degree + 1;
    Matrix normal(m, Vector(m, 0.0));
    Vector rhs(m, 0.0);

    for (std::size_t row = 0; row < m; ++row) {
        for (std::size_t col = 0; col < m; ++col) {
            for (double x : xs) {
                normal[row][col] += std::pow(x, static_cast<double>(row + col));
            }
        }
        for (std::size_t i = 0; i < xs.size(); ++i) {
            rhs[row] += ys[i] * std::pow(xs[i], static_cast<double>(row));
        }
    }

    return gaussian_elimination(normal, rhs).solution;
}

double evaluate_polynomial(
    const Vector& coefficients,
    double x)
{
    double result = 0.0;
    for (std::size_t i = coefficients.size(); i-- > 0;) {
        result = result * x + coefficients[i];
    }
    return result;
}

Vector chebyshev_nodes(double a, double b, std::size_t count)
{
    if (count == 0) {
        throw std::invalid_argument("count must be positive");
    }
    Vector nodes(count, 0.0);
    for (std::size_t k = 0; k < count; ++k) {
        const double angle = (2.0 * static_cast<double>(k) + 1.0) * kPi / (2.0 * static_cast<double>(count));
        nodes[k] = 0.5 * (a + b) + 0.5 * (b - a) * std::cos(angle);
    }
    std::sort(nodes.begin(), nodes.end());
    return nodes;
}

double chebyshev_polynomial(int degree, double x)
{
    if (degree < 0) {
        throw std::invalid_argument("degree must be non-negative");
    }
    if (degree == 0) {
        return 1.0;
    }
    if (degree == 1) {
        return x;
    }
    double t0 = 1.0;
    double t1 = x;
    for (int n = 1; n < degree; ++n) {
        const double t2 = 2.0 * x * t1 - t0;
        t0 = t1;
        t1 = t2;
    }
    return t1;
}

ComplexVector fft(const ComplexVector& values, bool inverse)
{
    const std::size_t n = values.size();
    if (n == 0 || (n & (n - 1)) != 0) {
        throw std::invalid_argument("FFT input size must be a positive power of two");
    }

    ComplexVector a = values;
    for (std::size_t i = 1, j = 0; i < n; ++i) {
        std::size_t bit = n >> 1;
        for (; j & bit; bit >>= 1) {
            j ^= bit;
        }
        j ^= bit;
        if (i < j) {
            std::swap(a[i], a[j]);
        }
    }

    for (std::size_t length = 2; length <= n; length <<= 1) {
        const double angle = 2.0 * kPi / static_cast<double>(length) * (inverse ? 1.0 : -1.0);
        const Complex root(std::cos(angle), std::sin(angle));
        for (std::size_t i = 0; i < n; i += length) {
            Complex w = 1.0;
            for (std::size_t j = 0; j < length / 2; ++j) {
                const Complex u = a[i + j];
                const Complex v = a[i + j + length / 2] * w;
                a[i + j] = u + v;
                a[i + j + length / 2] = u - v;
                w *= root;
            }
        }
    }

    if (inverse) {
        for (Complex& value : a) {
            value /= static_cast<double>(n);
        }
    }
    return a;
}

} // namespace numerical
