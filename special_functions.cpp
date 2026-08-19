#include "special_functions.hpp"

#include <algorithm>
#include <cmath>
#include <cstddef>

namespace numerical {
namespace {

constexpr double kSeriesTolerance = 1e-16;
constexpr std::size_t kMaxSeriesTerms = 200;

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

} // namespace numerical
