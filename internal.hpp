#pragma once

#include "errors.hpp"
#include "types.hpp"

#include <cmath>
#include <stdexcept>

namespace numerical::detail {

constexpr double kTiny = 1e-14;

inline void require_same_size_nonempty(const Vector& xs, const Vector& ys)
{
    if (xs.empty() || xs.size() != ys.size()) {
        throw std::invalid_argument("vectors must be non-empty and have the same size");
    }
}

inline void require_square(const Matrix& a)
{
    if (a.empty()) {
        throw std::invalid_argument("matrix must be non-empty");
    }
    for (const auto& row : a) {
        if (row.size() != a.size()) {
            throw std::invalid_argument("matrix must be square");
        }
    }
}

inline double l2_distance(const Vector& a, const Vector& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("vector dimensions are incompatible");
    }
    Vector difference(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        difference[i] = a[i] - b[i];
    }
    return l2_norm(difference);
}

inline Vector multiply(const Matrix& a, const Vector& x)
{
    Vector result(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        if (a[i].size() != x.size()) {
            throw std::invalid_argument("matrix and vector dimensions are incompatible");
        }
        for (std::size_t j = 0; j < x.size(); ++j) {
            result[i] += a[i][j] * x[j];
        }
    }
    return result;
}

inline Vector subtract(const Vector& a, const Vector& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("vector dimensions are incompatible");
    }
    Vector result(a.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        result[i] = a[i] - b[i];
    }
    return result;
}

inline double dot(const Vector& a, const Vector& b)
{
    if (a.size() != b.size()) {
        throw std::invalid_argument("vector dimensions are incompatible");
    }
    double sum = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) {
        sum += a[i] * b[i];
    }
    return sum;
}

} // namespace numerical::detail
