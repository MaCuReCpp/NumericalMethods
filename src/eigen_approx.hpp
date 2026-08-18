#pragma once

#include "types.hpp"

namespace numerical {

EigenResult power_method(
    const Matrix& a,
    Vector initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 1000);

QRDecomposition qr_decomposition(const Matrix& a);

Vector qr_eigenvalues(
    Matrix a,
    Tolerance tolerance = {},
    std::size_t max_iterations = 1000);

SVDResult svd(
    const Matrix& a,
    Tolerance tolerance = {},
    std::size_t max_iterations = 1000);

} // namespace numerical
