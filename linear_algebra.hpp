#pragma once

#include "types.hpp"

namespace numerical {

LinearSolveResult gaussian_elimination(Matrix a, Vector b);

LUDecomposition lu_decomposition(Matrix a);
Vector solve_lu(const LUDecomposition& lu, const Vector& b);
double determinant(Matrix a);
Matrix inverse(Matrix a);

LinearSolveResult jacobi(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 500);

LinearSolveResult gauss_seidel(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 500);

LinearSolveResult sor(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    double omega,
    Tolerance tolerance = {},
    std::size_t max_iterations = 500);

LinearSolveResult conjugate_gradient(
    const Matrix& a,
    const Vector& b,
    Vector initial_guess,
    Tolerance tolerance = {},
    std::size_t max_iterations = 1000);

LinearSolveResult iterative_refinement(
    const Matrix& a,
    const Vector& b,
    Vector initial_solution,
    Tolerance tolerance = {},
    std::size_t max_iterations = 10);

} // namespace numerical
