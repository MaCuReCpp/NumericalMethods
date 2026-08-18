#include "quantum.hpp"

#include "differentiation_integration.hpp"
#include "linear_algebra.hpp"

#include <cmath>
#include <stdexcept>

namespace numerical::quantum {
namespace {

constexpr double kTiny = 1e-14;
constexpr double kPi = 3.141592653589793238462643383279502884;
const Complex kI{0.0, 1.0};

void require_same_size(const ComplexVector& a, const ComplexVector& b)
{
    if (a.empty() || a.size() != b.size()) {
        throw std::invalid_argument("complex vectors must be non-empty and have the same size");
    }
}

void require_rectangular(const ComplexMatrix& matrix)
{
    if (matrix.empty() || matrix.front().empty()) {
        throw std::invalid_argument("matrix must be non-empty");
    }
    const std::size_t cols = matrix.front().size();
    for (const auto& row : matrix) {
        if (row.size() != cols) {
            throw std::invalid_argument("matrix must be rectangular");
        }
    }
}

void require_square(const ComplexMatrix& matrix)
{
    require_rectangular(matrix);
    if (matrix.size() != matrix.front().size()) {
        throw std::invalid_argument("matrix must be square");
    }
}

double grid_spacing(const Vector& x)
{
    if (x.size() < 2) {
        throw std::invalid_argument("grid must contain at least two points");
    }
    const double dx = x[1] - x[0];
    if (dx <= 0.0) {
        throw std::invalid_argument("grid must be strictly increasing");
    }
    for (std::size_t i = 1; i + 1 < x.size(); ++i) {
        if (std::abs((x[i + 1] - x[i]) - dx) > 1e-10 * std::max(1.0, std::abs(dx))) {
            throw std::invalid_argument("grid must be uniformly spaced");
        }
    }
    return dx;
}

double integrate_density(const Vector& x, const ComplexVector& values)
{
    if (x.size() != values.size() || x.size() < 2) {
        throw std::invalid_argument("wavefunction grid dimensions are incompatible");
    }
    double sum = 0.0;
    for (std::size_t i = 0; i + 1 < x.size(); ++i) {
        sum += 0.5 * (std::norm(values[i]) + std::norm(values[i + 1])) * (x[i + 1] - x[i]);
    }
    return sum;
}

} // namespace

Complex inner_product(const ComplexVector& bra, const ComplexVector& ket)
{
    require_same_size(bra, ket);
    Complex sum = 0.0;
    for (std::size_t i = 0; i < bra.size(); ++i) {
        sum += std::conj(bra[i]) * ket[i];
    }
    return sum;
}

double norm(const ComplexVector& state)
{
    if (state.empty()) {
        throw std::invalid_argument("state must be non-empty");
    }
    return std::sqrt(std::real(inner_product(state, state)));
}

ComplexVector normalize(const ComplexVector& state)
{
    const double n = norm(state);
    if (n < kTiny) {
        throw std::invalid_argument("cannot normalize a zero state");
    }
    ComplexVector result = state;
    for (auto& value : result) {
        value /= n;
    }
    return result;
}

double probability(const ComplexVector& state, std::size_t index)
{
    const ComplexVector normalized = normalize(state);
    if (index >= normalized.size()) {
        throw std::out_of_range("basis index is out of range");
    }
    return std::norm(normalized[index]);
}

Vector probability_distribution(const ComplexVector& state)
{
    const ComplexVector normalized = normalize(state);
    Vector probabilities(normalized.size(), 0.0);
    for (std::size_t i = 0; i < normalized.size(); ++i) {
        probabilities[i] = std::norm(normalized[i]);
    }
    return probabilities;
}

ComplexMatrix identity(std::size_t dimension)
{
    ComplexMatrix result(dimension, ComplexVector(dimension, 0.0));
    for (std::size_t i = 0; i < dimension; ++i) {
        result[i][i] = 1.0;
    }
    return result;
}

ComplexMatrix dagger(const ComplexMatrix& matrix)
{
    require_rectangular(matrix);
    ComplexMatrix result(matrix.front().size(), ComplexVector(matrix.size(), 0.0));
    for (std::size_t i = 0; i < matrix.size(); ++i) {
        for (std::size_t j = 0; j < matrix.front().size(); ++j) {
            result[j][i] = std::conj(matrix[i][j]);
        }
    }
    return result;
}

ComplexMatrix add(const ComplexMatrix& a, const ComplexMatrix& b)
{
    require_rectangular(a);
    require_rectangular(b);
    if (a.size() != b.size() || a.front().size() != b.front().size()) {
        throw std::invalid_argument("matrix dimensions are incompatible");
    }
    ComplexMatrix result = a;
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < a.front().size(); ++j) {
            result[i][j] += b[i][j];
        }
    }
    return result;
}

ComplexMatrix subtract(const ComplexMatrix& a, const ComplexMatrix& b)
{
    return add(a, scalar_multiply(-1.0, b));
}

ComplexMatrix multiply(const ComplexMatrix& a, const ComplexMatrix& b)
{
    require_rectangular(a);
    require_rectangular(b);
    if (a.front().size() != b.size()) {
        throw std::invalid_argument("matrix dimensions are incompatible");
    }
    ComplexMatrix result(a.size(), ComplexVector(b.front().size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t k = 0; k < b.size(); ++k) {
            for (std::size_t j = 0; j < b.front().size(); ++j) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
    return result;
}

ComplexVector apply(const ComplexMatrix& operator_matrix, const ComplexVector& state)
{
    require_rectangular(operator_matrix);
    if (operator_matrix.front().size() != state.size()) {
        throw std::invalid_argument("operator and state dimensions are incompatible");
    }
    ComplexVector result(operator_matrix.size(), 0.0);
    for (std::size_t i = 0; i < operator_matrix.size(); ++i) {
        for (std::size_t j = 0; j < state.size(); ++j) {
            result[i] += operator_matrix[i][j] * state[j];
        }
    }
    return result;
}

ComplexMatrix scalar_multiply(Complex scalar, const ComplexMatrix& matrix)
{
    require_rectangular(matrix);
    ComplexMatrix result = matrix;
    for (auto& row : result) {
        for (auto& value : row) {
            value *= scalar;
        }
    }
    return result;
}

ComplexMatrix commutator(const ComplexMatrix& a, const ComplexMatrix& b)
{
    return subtract(multiply(a, b), multiply(b, a));
}

ComplexMatrix tensor_product(const ComplexMatrix& a, const ComplexMatrix& b)
{
    require_rectangular(a);
    require_rectangular(b);
    ComplexMatrix result(a.size() * b.size(), ComplexVector(a.front().size() * b.front().size(), 0.0));
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < a.front().size(); ++j) {
            for (std::size_t k = 0; k < b.size(); ++k) {
                for (std::size_t l = 0; l < b.front().size(); ++l) {
                    result[i * b.size() + k][j * b.front().size() + l] = a[i][j] * b[k][l];
                }
            }
        }
    }
    return result;
}

ComplexVector tensor_product(const ComplexVector& a, const ComplexVector& b)
{
    if (a.empty() || b.empty()) {
        throw std::invalid_argument("states must be non-empty");
    }
    ComplexVector result(a.size() * b.size(), 0.0);
    for (std::size_t i = 0; i < a.size(); ++i) {
        for (std::size_t j = 0; j < b.size(); ++j) {
            result[i * b.size() + j] = a[i] * b[j];
        }
    }
    return result;
}

Complex expectation_value(const ComplexMatrix& operator_matrix, const ComplexVector& state)
{
    require_square(operator_matrix);
    const ComplexVector normalized = normalize(state);
    return inner_product(normalized, apply(operator_matrix, normalized));
}

double variance(const ComplexMatrix& operator_matrix, const ComplexVector& state)
{
    const Complex mean = expectation_value(operator_matrix, state);
    const ComplexMatrix squared = multiply(operator_matrix, operator_matrix);
    const Complex second = expectation_value(squared, state);
    return std::max(0.0, std::real(second - mean * mean));
}

double uncertainty(const ComplexMatrix& operator_matrix, const ComplexVector& state)
{
    return std::sqrt(variance(operator_matrix, state));
}

Complex transition_amplitude(const ComplexVector& from, const ComplexVector& to)
{
    return inner_product(normalize(to), normalize(from));
}

ComplexVector project_onto(const ComplexVector& state, const ComplexVector& basis_state)
{
    const ComplexVector basis = normalize(basis_state);
    const Complex coefficient = inner_product(basis, state);
    ComplexVector projection = basis;
    for (auto& value : projection) {
        value *= coefficient;
    }
    return projection;
}

ComplexMatrix pauli_x()
{
    return {{0.0, 1.0}, {1.0, 0.0}};
}

ComplexMatrix pauli_y()
{
    return {{0.0, -kI}, {kI, 0.0}};
}

ComplexMatrix pauli_z()
{
    return {{1.0, 0.0}, {0.0, -1.0}};
}

ComplexMatrix spin_x(double hbar)
{
    return scalar_multiply(0.5 * hbar, pauli_x());
}

ComplexMatrix spin_y(double hbar)
{
    return scalar_multiply(0.5 * hbar, pauli_y());
}

ComplexMatrix spin_z(double hbar)
{
    return scalar_multiply(0.5 * hbar, pauli_z());
}

ComplexMatrix spin_plus(double hbar)
{
    return {{0.0, hbar}, {0.0, 0.0}};
}

ComplexMatrix spin_minus(double hbar)
{
    return {{0.0, 0.0}, {hbar, 0.0}};
}

double infinite_square_well_energy(int n, double mass, double length, double hbar)
{
    if (n <= 0 || mass <= 0.0 || length <= 0.0 || hbar <= 0.0) {
        throw std::invalid_argument("invalid infinite square well parameters");
    }
    return (static_cast<double>(n * n) * kPi * kPi * hbar * hbar) / (2.0 * mass * length * length);
}

double infinite_square_well_wavefunction(int n, double x, double length)
{
    if (n <= 0 || length <= 0.0) {
        throw std::invalid_argument("invalid infinite square well parameters");
    }
    if (x < 0.0 || x > length) {
        return 0.0;
    }
    return std::sqrt(2.0 / length) * std::sin(static_cast<double>(n) * kPi * x / length);
}

double harmonic_oscillator_energy(int n, double omega, double hbar)
{
    if (n < 0 || omega <= 0.0 || hbar <= 0.0) {
        throw std::invalid_argument("invalid harmonic oscillator parameters");
    }
    return hbar * omega * (static_cast<double>(n) + 0.5);
}

ComplexMatrix harmonic_oscillator_annihilation(std::size_t dimension)
{
    ComplexMatrix result(dimension, ComplexVector(dimension, 0.0));
    for (std::size_t n = 1; n < dimension; ++n) {
        result[n - 1][n] = std::sqrt(static_cast<double>(n));
    }
    return result;
}

ComplexMatrix harmonic_oscillator_creation(std::size_t dimension)
{
    return dagger(harmonic_oscillator_annihilation(dimension));
}

ComplexMatrix harmonic_oscillator_number(std::size_t dimension)
{
    return multiply(harmonic_oscillator_creation(dimension), harmonic_oscillator_annihilation(dimension));
}

Wavefunction1D normalize_wavefunction(const Wavefunction1D& wavefunction)
{
    const double integral = integrate_density(wavefunction.x, wavefunction.psi);
    if (integral < kTiny) {
        throw std::invalid_argument("cannot normalize a zero wavefunction");
    }
    Wavefunction1D result = wavefunction;
    const double factor = std::sqrt(integral);
    for (auto& value : result.psi) {
        value /= factor;
    }
    return result;
}

double position_expectation(const Wavefunction1D& wavefunction)
{
    const auto normalized = normalize_wavefunction(wavefunction);
    double sum = 0.0;
    for (std::size_t i = 0; i + 1 < normalized.x.size(); ++i) {
        const double left = normalized.x[i] * std::norm(normalized.psi[i]);
        const double right = normalized.x[i + 1] * std::norm(normalized.psi[i + 1]);
        sum += 0.5 * (left + right) * (normalized.x[i + 1] - normalized.x[i]);
    }
    return sum;
}

double momentum_expectation(const Wavefunction1D& wavefunction, double hbar)
{
    if (hbar <= 0.0) {
        throw std::invalid_argument("hbar must be positive");
    }
    const auto normalized = normalize_wavefunction(wavefunction);
    const double dx = grid_spacing(normalized.x);
    Complex integral = 0.0;
    for (std::size_t i = 1; i + 1 < normalized.x.size(); ++i) {
        const Complex derivative = (normalized.psi[i + 1] - normalized.psi[i - 1]) / (2.0 * dx);
        integral += std::conj(normalized.psi[i]) * (-kI * hbar) * derivative * dx;
    }
    return std::real(integral);
}

DiscreteHamiltonian1D finite_difference_hamiltonian(
    const Vector& x,
    const Vector& potential,
    double mass,
    double hbar)
{
    if (x.size() != potential.size() || x.size() < 3 || mass <= 0.0 || hbar <= 0.0) {
        throw std::invalid_argument("invalid Hamiltonian grid or physical parameters");
    }
    const double dx = grid_spacing(x);
    const std::size_t n = x.size();
    Matrix hamiltonian(n, Vector(n, 0.0));
    const double kinetic_diag = hbar * hbar / (mass * dx * dx);
    const double kinetic_offdiag = -hbar * hbar / (2.0 * mass * dx * dx);
    for (std::size_t i = 0; i < n; ++i) {
        hamiltonian[i][i] = kinetic_diag + potential[i];
        if (i > 0) {
            hamiltonian[i][i - 1] = kinetic_offdiag;
        }
        if (i + 1 < n) {
            hamiltonian[i][i + 1] = kinetic_offdiag;
        }
    }
    return {x, hamiltonian};
}

} // namespace numerical::quantum
