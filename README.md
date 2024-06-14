# Quantum Numerical Methods Project

## Overview

This project focuses on the implementation and analysis of various quantum numerical methods. It encompasses exact diagonalization, polynomial interpolation, spline interpolation, and the numerical analysis of key quantum Hamiltonians such as the Bose-Hubbard Hamiltonian and the XX-Z Hamiltonian. Additionally, it includes a numerical approximation for the density matrix formulation used in quantum renormalization groups and methods to solve the Gross-Pitaevskii equation.

## Table of Contents

1. [Project Structure](#project-structure)
2. [Exact Diagonalization](#exact-diagonalization)
3. [Polynomial Interpolation](#polynomial-interpolation)
4. [Spline Interpolation](#spline-interpolation)
5. [Numerical Analysis of Hamiltonians](#numerical-analysis-of-hamiltonians)
   - [Bose-Hubbard Hamiltonian](#bose-hubbard-hamiltonian)
   - [XX-Z Hamiltonian](#xx-z-hamiltonian)
6. [Density Matrix Renormalization Group (DMRG)](#density-matrix-renormalization-group-dmrg)
7. [Solving the Gross-Pitaevskii Equation](#solving-the-gross-pitaevskii-equation)
8. [Dependencies](#dependencies)
9. [Installation](#installation)
10. [Usage](#usage)
11. [Contributing](#contributing)

## Project Structure

The project directory is organized as follows:

quantum-numerical-methods/
├── README.md
├── HW_1/
│ ├──problem_1
│ │ ├── least_squares.py
│ │ ├── polynomial_interpolation.py
│ │ ├── spline_interpolation.py
│ ├──problem_2
│ │ ├── xxz_hamiltonian.py
│ ├──problem_3
│ ├── bose_hubbard.py
├── HW_2/
│ ├── dmrg.py
│ └── gross_pitaevskii.py
└── requirements.txt

## Exact Diagonalization

Exact diagonalization is a method used to find the eigenvalues and eigenvectors of Hamiltonians. This method is computationally intensive but provides precise solutions for small systems. The `exact_diagonalization.py` module includes functions to:

- Construct the Hamiltonian matrix.
- Perform exact diagonalization using linear algebra libraries.
- Analyze the spectrum and eigenstates.

## Polynomial Interpolation

Polynomial interpolation is used to estimate values between known data points. The `polynomial_interpolation.py` module includes functions to:

- Construct interpolating polynomials.
- Evaluate the polynomial at given points.
- Analyze interpolation errors.

## Spline Interpolation

Spline interpolation uses piecewise polynomials, typically cubic splines, for a smoother interpolation than polynomial interpolation. The `spline_interpolation.py` module includes functions to:

- Construct spline functions.
- Evaluate splines at given points.
- Analyze spline fitting and interpolation errors.

## Numerical Analysis of Hamiltonians

### Bose-Hubbard Hamiltonian

The Bose-Hubbard model describes interacting bosons on a lattice. The `bose_hubbard.py` module includes functions to:

- Construct the Bose-Hubbard Hamiltonian matrix.
- Perform numerical diagonalization.
- Analyze the ground state and excitation spectrum.

### XX-Z Hamiltonian

The XX-Z model describes spin-1/2 particles on a lattice with anisotropic interactions. The `xxz_hamiltonian.py` module includes functions to:

- Construct the XX-Z Hamiltonian matrix.
- Perform numerical diagonalization.
- Analyze the magnetic properties and phase transitions.

## Density Matrix Renormalization Group (DMRG)

DMRG is a numerical variational technique designed to study low-dimensional quantum systems. The `dmrg.py` module includes functions to:

- Initialize the density matrix renormalization group algorithm.
- Perform iterative optimization to approximate the ground state.
- Analyze the convergence and accuracy of the results.

## Solving the Gross-Pitaevskii Equation

The Gross-Pitaevskii equation describes the wave function of a Bose-Einstein condensate. The `gross_pitaevskii.py` module includes functions to:

- Discretize the Gross-Pitaevskii equation.
- Implement numerical solvers (e.g., finite difference, Runge-Kutta).
- Analyze the dynamics and stationary solutions of the condensate.

## Dependencies

The project requires the following Python libraries:

- numpy
- scipy
- matplotlib

## Installation

To install the dependencies, run:

pip install -r requirements.txt

vbnet

## Usage

Each module can be executed independently. Example usage is provided within each module as well as in the `tests` directory. To run the tests, execute:

pytest tests/

css

## Contributing

Contributions are welcome! Please fork the repository and create a pull request with your changes. Ensure that your code passes all tests and adheres to the project's coding standards.
