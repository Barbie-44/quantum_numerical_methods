import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.linalg import svd
import matplotlib.pyplot as plt
import numpy as np
import scipy.linalg as la


class CaldeiraLeggett:
    """
    GOAL:
        Get eigenvalues and eigenstates of the Caldeira-Legget Hamiltonian via Exact Diagonalization(ED).
    REMARKS:
        .By simplicity all masses and Planck constant are equal to 1.
        .The number of bath oscillators is finite.
    """

    def __init__(self, dimension_sys, N_reservoirs, dimension_reservoir):
        self.omega_sys = 1.0
        self.dim_sys = dimension_sys

        # Parameters of reservoir
        self.N_bath = N_reservoirs
        self.omega_bath = np.linspace(
            0.5, 2.0, self.N_bath
        )  # frequencies of the bath oscillators
        self.dim_bath = dimension_reservoir

        # Coupling constants for system-bath interaction
        self.c_bath = np.ones(self.N_bath) * 0.01  # coupling constants

        # Total Hilbert space dimension (system and bath)
        self.total_dim = self.dim_sys * (self.dim_bath**self.N_bath)

    # Function to create the creation and annihilation operators in the Fock basis
    def creation_operator(self, dim):
        """Returns the creation operator for a harmonic oscillator with dimension dim"""
        return sp.diags(np.sqrt(np.arange(1, dim)), 1, format="csr")

    def annihilation_operator(self, dim):
        """Returns the annihilation operator for a harmonic oscillator with dimension dim"""
        return sp.diags(np.sqrt(np.arange(1, dim)), -1, format="csr")

    # Position operator in terms of creation and annihilation operators: q = (a + a†) / sqrt(2)
    def position_operator(self, dim):
        a = self.annihilation_operator(dim)
        a_dag = self.creation_operator(dim)
        return (a + a_dag) / np.sqrt(2)

    # Create the system Hamiltonian (harmonic oscillator in Fock basis)
    def system_hamiltonian(self, dim_sys, omega_sys):
        diag_elements = [(n + 0.5) * omega_sys for n in range(dim_sys)]
        return sp.diags(diag_elements, format="csr")

    # Create the bath Hamiltonian for a single bath oscillator
    def single_bath_hamiltonian(self, dim_bath, omega):
        diag_elements = [(n + 0.5) * omega for n in range(dim_bath)]
        return sp.diags(diag_elements, format="csr")

    # Full bath Hamiltonian (sum over all bath oscillators)
    def bath_hamiltonian(self, dim_bath, omega_bath):
        H_bath_total = sp.csr_matrix(
            (dim_bath**self.N_bath, dim_bath**self.N_bath)
        )  # Initialize the full bath Hamiltonian
        for b, omega in enumerate(omega_bath):
            H_b = self.single_bath_hamiltonian(dim_bath, omega)
            # Kronecker product to apply H_b to the correct bath oscillator
            identity_before = sp.eye(dim_bath**b)
            identity_after = sp.eye(dim_bath ** (self.N_bath - b - 1))
            H_bath_total += sp.kron(
                sp.kron(identity_before, H_b), identity_after
            )
        return H_bath_total

    # Create the system-bath interaction term (linear coupling: q_sys * q_bath)
    def interaction_hamiltonian(self, dim_sys, dim_bath, c_bath):
        q_sys = self.position_operator(
            dim_sys
        )  # Position operator for the system
        H_int = sp.lil_matrix(
            (self.total_dim, self.total_dim)
        )  # Use LIL for efficient element-wise manipulation
        for b, c in enumerate(c_bath):
            q_bath = self.position_operator(
                dim_bath
            )  # Position operator for the bath
            identity_before = sp.eye(dim_bath**b)
            identity_after = sp.eye(dim_bath ** (self.N_bath - b - 1))
            # q_sys interacts with each bath oscillator q_bath
            H_int += c * sp.kron(
                q_sys,
                sp.kron(identity_before, sp.kron(q_bath, identity_after)),
            )
        return H_int.tocsr()  # Convert to CSR format for diagonalization

    # Total Hamiltonian: system + bath + interaction
    def total_hamiltonian(self, H_sys, H_bath, H_int):
        H_total = (
            sp.kron(H_sys, sp.eye(self.dim_bath**self.N_bath))
            + sp.kron(sp.eye(self.dim_sys), H_bath)
            + H_int
        )
        return H_total

    def get_total_hamiltonian(self):
        # Construct the system Hamiltonian
        H_sys = self.system_hamiltonian(self.dim_sys, self.omega_sys)
        print("H SYS: ")
        print(H_sys)

        # Construct the bath Hamiltonian
        H_bath = self.bath_hamiltonian(self.dim_bath, self.omega_bath)
        print("H BATH: ")
        print(H_bath)

        # Construct the interaction Hamiltonian (linear coupling)
        H_int = self.interaction_hamiltonian(
            self.dim_sys, self.dim_bath, self.c_bath
        )
        print("H INT")
        print(H_int)

        # Construct the total Hamiltonian
        H_total = self.total_hamiltonian(H_sys, H_bath, H_int)
        print("H TOTAL")
        print(H_total)
        return H_total

    # Diagonalize the Hamiltonian (sparse)
    def diagonalize_sparse(self, H_total, k=1):
        eigvals, eigvecs = spla.eigsh(
            H_total, k=1, which="SM"
        )  # Smallest k eigenvalues/eigenvectors
        return eigvals, eigvecs

    def calculate_trace_norm_transpose_matrix(self, rho):
        # H_total = self.get_total_hamiltonian()
        # eigvals, eigvecs = self.diagonalize_sparse(H_total)
        # print("Eigenvalues of the total Hamiltonian:\n", eigvals)
        # print("Eigenstates of the total Hamiltonian:\n", eigvecs)
        # ground_state = eigvecs[:, 0]
        # rho = np.outer(ground_state, ground_state.conj())
        # print("RHO: ", rho)
        rho_transposed = np.transpose(rho)
        singular_values = svd(rho_transposed, compute_uv=False)
        trace_norm = np.sum(np.abs(singular_values))
        print("Trace Norm of the partially transposed matrix:", trace_norm)
        return trace_norm

    def calculate_logarithmic_negativity(self, rho):
        trace_rho = self.calculate_trace_norm_transpose_matrix(rho)
        E_N = np.log(trace_rho)
        print("EN: ", E_N)
        return E_N


dimension_sys = 2
N_reservoirs = 10
dimension_reservoir = 2
result = CaldeiraLeggett(dimension_sys, 2, dimension_reservoir)
# result.calculate_logarithmic_negativity()
H_total = result.get_total_hamiltonian()


# Constants
hbar = 1.0545718e-34  # Reduced Planck's constant, in Js
t0 = 0  # Initial time
t = np.linspace(t0, 60, 200)  # Final time, example in seconds
eigenvalues, eigenvectors = result.diagonalize_sparse(H_total)
ground_state = eigenvectors[:, 0]
rho0 = np.outer(ground_state, ground_state.conj())  # Pure state density matrix


def time_evolution_operator(eigenvalues, t, t0, hbar):
    # U = diag(exp(-i * E_i * (t - t0) / hbar))
    phase_factors = np.exp(-1j * eigenvalues * (t - t0) / hbar)
    return np.diag(phase_factors)


negativities = []
for t_end in t:
    U = time_evolution_operator(eigenvalues, t_end, t0, hbar)

    rho0_eigen = np.dot(eigenvectors.T.conj(), np.dot(rho0, eigenvectors))

    # Evolve the density matrix in the eigenbasis: rho(t)_eigen = U * rho0_eigen * U^dagger
    rho_t_eigen = np.dot(U, np.dot(rho0_eigen, U.conj().T))

    # Transform the density matrix back to the original basis: rho(t) = V * rho(t)_eigen * V^dagger
    rho_t = np.dot(eigenvectors, np.dot(rho_t_eigen, eigenvectors.T.conj()))

    # Output the time-evolved density matrix
    print("Density matrix at time t:")
    print(rho_t)
    E_N = result.calculate_logarithmic_negativity(rho_t)
    negativities.append(E_N)


# results = []
# for i in range(2, N_reservoirs):
#     print("*" * 100)
#     print("ITERATION:{i} " * 10)
#     result = CaldeiraLeggett(dimension_sys, i, dimension_reservoir)
#     negativity = result.calculate_logarithmic_negativity()
#     print("NEGATIVITY: ", negativity)
#     if negativity < 0:
#         negativity = 0
#     results.append(negativity)
# N = [N for N in range(2, N_reservoirs)]
plt.plot(t, negativities, "-")
plt.xlabel("t")
plt.ylabel("EN")
plt.title("Logaritmic negativity")
plt.legend()
plt.show()
