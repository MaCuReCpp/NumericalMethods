#include "interpolation.hpp"
#include "internal.hpp"

#include <algorithm>
#include <cmath>

namespace numerical {
using detail::require_same_size_nonempty;

std::size_t spline_interval(const CubicSpline& spline, double x)
{
    if (spline.nodes.size() < 2) {
        throw std::invalid_argument("spline must contain at least two nodes");
    }
    if (x <= spline.nodes.front()) {
        return 0;
    }
    if (x >= spline.nodes.back()) {
        return spline.nodes.size() - 2;
    }
    return static_cast<std::size_t>(
        std::upper_bound(spline.nodes.begin(), spline.nodes.end(), x) - spline.nodes.begin() - 1);
}


double neville_interpolate(
    const Vector& xs,
    const Vector& ys,
    double x)
{
    require_same_size_nonempty(xs, ys);
    Vector q = ys;
    for (std::size_t i = 1; i < xs.size(); ++i) {
        for (std::size_t j = 0; j < xs.size() - i; ++j) {
            const double denominator = xs[j] - xs[j + i];
            if (std::abs(denominator) < detail::kTiny) {
                throw std::invalid_argument("interpolation nodes must be distinct");
            }
            q[j] = ((x - xs[j + i]) * q[j] + (xs[j] - x) * q[j + 1]) / denominator;
        }
    }
    return q.front();
}

double hermite_interpolate(
    const Vector& xs,
    const Vector& ys,
    const Vector& derivatives,
    double x)
{
    require_same_size_nonempty(xs, ys);
    if (derivatives.size() != xs.size()) {
        throw std::invalid_argument("derivative vector dimension is incompatible");
    }

    const std::size_t n = xs.size();
    Vector z(2 * n, 0.0);
    Matrix q(2 * n, Vector(2 * n, 0.0));

    for (std::size_t i = 0; i < n; ++i) {
        z[2 * i] = xs[i];
        z[2 * i + 1] = xs[i];
        q[2 * i][0] = ys[i];
        q[2 * i + 1][0] = ys[i];
        q[2 * i + 1][1] = derivatives[i];
        if (i == 0) {
            q[2 * i][1] = derivatives[i];
        } else {
            q[2 * i][1] = (q[2 * i][0] - q[2 * i - 1][0]) / (z[2 * i] - z[2 * i - 1]);
        }
    }

    for (std::size_t i = 2; i < 2 * n; ++i) {
        for (std::size_t j = 2; j <= i; ++j) {
            q[i][j] = (q[i][j - 1] - q[i - 1][j - 1]) / (z[i] - z[i - j]);
        }
    }

    double result = q[2 * n - 1][2 * n - 1];
    for (std::size_t i = 2 * n - 1; i-- > 0;) {
        result = result * (x - z[i]) + q[i][i];
    }
    return result;
}

CubicSpline natural_cubic_spline(
    const Vector& xs,
    const Vector& ys)
{
    require_same_size_nonempty(xs, ys);
    if (xs.size() < 2) {
        throw std::invalid_argument("at least two nodes are required");
    }

    const std::size_t n = xs.size() - 1;
    Vector h(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        h[i] = xs[i + 1] - xs[i];
        if (h[i] <= 0.0) {
            throw std::invalid_argument("spline nodes must be strictly increasing");
        }
    }

    Vector alpha(n, 0.0);
    for (std::size_t i = 1; i < n; ++i) {
        alpha[i] = 3.0 * (ys[i + 1] - ys[i]) / h[i] - 3.0 * (ys[i] - ys[i - 1]) / h[i - 1];
    }

    Vector l(n + 1, 1.0);
    Vector mu(n + 1, 0.0);
    Vector z(n + 1, 0.0);
    for (std::size_t i = 1; i < n; ++i) {
        l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }

    CubicSpline spline;
    spline.nodes = xs;
    spline.a = ys;
    spline.b.assign(n, 0.0);
    spline.c.assign(n + 1, 0.0);
    spline.d.assign(n, 0.0);

    for (std::size_t j = n; j-- > 0;) {
        spline.c[j] = z[j] - mu[j] * spline.c[j + 1];
        spline.b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (spline.c[j + 1] + 2.0 * spline.c[j]) / 3.0;
        spline.d[j] = (spline.c[j + 1] - spline.c[j]) / (3.0 * h[j]);
    }

    return spline;
}

CubicSpline clamped_cubic_spline(
    const Vector& xs,
    const Vector& ys,
    double left_derivative,
    double right_derivative)
{
    require_same_size_nonempty(xs, ys);
    if (xs.size() < 2) {
        throw std::invalid_argument("at least two nodes are required");
    }

    const std::size_t n = xs.size() - 1;
    Vector h(n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
        h[i] = xs[i + 1] - xs[i];
        if (h[i] <= 0.0) {
            throw std::invalid_argument("spline nodes must be strictly increasing");
        }
    }

    Vector alpha(n + 1, 0.0);
    alpha[0] = 3.0 * (ys[1] - ys[0]) / h[0] - 3.0 * left_derivative;
    alpha[n] = 3.0 * right_derivative - 3.0 * (ys[n] - ys[n - 1]) / h[n - 1];
    for (std::size_t i = 1; i < n; ++i) {
        alpha[i] = 3.0 * (ys[i + 1] - ys[i]) / h[i] - 3.0 * (ys[i] - ys[i - 1]) / h[i - 1];
    }

    Vector l(n + 1, 0.0);
    Vector mu(n + 1, 0.0);
    Vector z(n + 1, 0.0);
    l[0] = 2.0 * h[0];
    mu[0] = 0.5;
    z[0] = alpha[0] / l[0];
    for (std::size_t i = 1; i < n; ++i) {
        l[i] = 2.0 * (xs[i + 1] - xs[i - 1]) - h[i - 1] * mu[i - 1];
        mu[i] = h[i] / l[i];
        z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
    }
    l[n] = h[n - 1] * (2.0 - mu[n - 1]);
    z[n] = (alpha[n] - h[n - 1] * z[n - 1]) / l[n];

    CubicSpline spline;
    spline.nodes = xs;
    spline.a = ys;
    spline.b.assign(n, 0.0);
    spline.c.assign(n + 1, 0.0);
    spline.d.assign(n, 0.0);
    spline.c[n] = z[n];

    for (std::size_t j = n; j-- > 0;) {
        spline.c[j] = z[j] - mu[j] * spline.c[j + 1];
        spline.b[j] = (ys[j + 1] - ys[j]) / h[j] - h[j] * (spline.c[j + 1] + 2.0 * spline.c[j]) / 3.0;
        spline.d[j] = (spline.c[j + 1] - spline.c[j]) / (3.0 * h[j]);
    }
    return spline;
}

double evaluate_spline(const CubicSpline& spline, double x)
{
    const std::size_t i = spline_interval(spline, x);
    const double dx = x - spline.nodes[i];
    return spline.a[i] + spline.b[i] * dx + spline.c[i] * dx * dx + spline.d[i] * dx * dx * dx;
}

} // namespace numerical
