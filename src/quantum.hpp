#pragma once

#include "types.hpp"

namespace numerical::quantum {

struct QuantumState {
    ComplexVector amplitudes;
};

struct Wavefunction1D {
    Vector x;
    ComplexVector psi;
};

struct DiscreteHamiltonian1D {
    Vector x;
    Matrix matrix;
};

Complex inner_product(const ComplexVector& bra, const ComplexVector& ket);
double norm(const ComplexVector& state);
ComplexVector normalize(const ComplexVector& state);
double probability(const ComplexVector& state, std::size_t index);
Vector probability_distribution(const ComplexVector& state);

ComplexMatrix identity(std::size_t dimension);
ComplexMatrix dagger(const ComplexMatrix& matrix);
ComplexMatrix add(const ComplexMatrix& a, const ComplexMatrix& b);
ComplexMatrix subtract(const ComplexMatrix& a, const ComplexMatrix& b);
ComplexMatrix multiply(const ComplexMatrix& a, const ComplexMatrix& b);
ComplexVector apply(const ComplexMatrix& operator_matrix, const ComplexVector& state);
ComplexMatrix scalar_multiply(Complex scalar, const ComplexMatrix& matrix);
ComplexMatrix commutator(const ComplexMatrix& a, const ComplexMatrix& b);
ComplexMatrix tensor_product(const ComplexMatrix& a, const ComplexMatrix& b);
ComplexVector tensor_product(const ComplexVector& a, const ComplexVector& b);

Complex expectation_value(const ComplexMatrix& operator_matrix, const ComplexVector& state);
double variance(const ComplexMatrix& operator_matrix, const ComplexVector& state);
double uncertainty(const ComplexMatrix& operator_matrix, const ComplexVector& state);
Complex transition_amplitude(const ComplexVector& from, const ComplexVector& to);
ComplexVector project_onto(const ComplexVector& state, const ComplexVector& basis_state);

ComplexMatrix pauli_x();
ComplexMatrix pauli_y();
ComplexMatrix pauli_z();
ComplexMatrix spin_x(double hbar = 1.0);
ComplexMatrix spin_y(double hbar = 1.0);
ComplexMatrix spin_z(double hbar = 1.0);
ComplexMatrix spin_plus(double hbar = 1.0);
ComplexMatrix spin_minus(double hbar = 1.0);

double infinite_square_well_energy(int n, double mass, double length, double hbar = 1.0);
double infinite_square_well_wavefunction(int n, double x, double length);
double harmonic_oscillator_energy(int n, double omega, double hbar = 1.0);
ComplexMatrix harmonic_oscillator_annihilation(std::size_t dimension);
ComplexMatrix harmonic_oscillator_creation(std::size_t dimension);
ComplexMatrix harmonic_oscillator_number(std::size_t dimension);

Wavefunction1D normalize_wavefunction(const Wavefunction1D& wavefunction);
double position_expectation(const Wavefunction1D& wavefunction);
double momentum_expectation(const Wavefunction1D& wavefunction, double hbar = 1.0);
DiscreteHamiltonian1D finite_difference_hamiltonian(
    const Vector& x,
    const Vector& potential,
    double mass,
    double hbar = 1.0);

} // namespace numerical::quantum
