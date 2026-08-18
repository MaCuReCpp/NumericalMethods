#pragma once

#include <complex>
#include <cstddef>
#include <functional>
#include <vector>

namespace numerical {

using ScalarFunction = std::function<double(double)>;
using OdeFunction = std::function<double(double, double)>;
using Vector = std::vector<double>;
using Matrix = std::vector<Vector>;
using Complex = std::complex<double>;
using ComplexVector = std::vector<Complex>;
using ComplexMatrix = std::vector<ComplexVector>;
using VectorFunction = std::function<Vector(const Vector&)>;
using JacobianFunction = std::function<Matrix(const Vector&)>;

struct Tolerance {
    double absolute = 1e-10;
    double relative = 1e-10;
};

struct IterationResult {
    double value = 0.0;
    double residual = 0.0;
    std::size_t iterations = 0;
    bool converged = false;
};

struct LinearSolveResult {
    Vector solution;
    double residual_norm = 0.0;
    std::size_t iterations = 0;
    bool converged = false;
};

struct OdePoint {
    double t = 0.0;
    double y = 0.0;
};

struct CubicSpline {
    Vector nodes;
    Vector a;
    Vector b;
    Vector c;
    Vector d;
};

struct LUDecomposition {
    Matrix lower;
    Matrix upper;
    std::vector<std::size_t> permutation;
    int row_swaps = 0;
};

struct EigenResult {
    double eigenvalue = 0.0;
    Vector eigenvector;
    double residual_norm = 0.0;
    std::size_t iterations = 0;
    bool converged = false;
};

struct Grid1D {
    Vector points;
    Vector values;
};

struct HeatEquationResult {
    Vector x;
    Vector t;
    Matrix values;
};

struct Grid2D {
    Vector x;
    Vector y;
    Matrix values;
};

struct WaveEquationResult {
    Vector x;
    Vector t;
    Matrix values;
};

struct QRDecomposition {
    Matrix q;
    Matrix r;
};

struct SVDResult {
    Matrix u;
    Vector singular_values;
    Matrix v_transpose;
};

enum class DifferenceMethod {
    Forward,
    Backward,
    Central
};

} // namespace numerical
