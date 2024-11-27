#include "Hamiltonnian.hpp"
#include <iostream>
#include <omp.h>
#include <complex>

// Constructor
Hamiltoniano::Hamiltoniano(int N, double J, double G)
    : N(N), J(J), G(G) {
    inicializar_matrices_pauli();
    hamiltoniano = std::vector<std::vector<std::complex<double>>>(
        1 << N, std::vector<std::complex<double>>(1 << N, {0.0, 0.0}));
}

// Inicialización de matrices de Pauli y la identidad
void Hamiltoniano::inicializar_matrices_pauli() {
    pauli_x = {{0, 1}, {1, 0}};
    pauli_z = {{1, 0}, {0, -1}};
    identidad = {{1, 0}, {0, 1}};
}

// Mostrar el Hamiltoniano
void Hamiltoniano::mostrar_Hamiltoniano() const {
    for (const auto& fila : hamiltoniano) {
        for (const auto& elem : fila) {
            std::cout << elem << " ";
        }
        std::cout << "\n";
    }
}

// Construcción del Hamiltoniano con paralelización
void Hamiltoniano::crear_Hamiltoniano() {
    // Parte relacionada con Pauli-Z
    #pragma omp parallel for
    for (int i = 0; i < N - 1; ++i) {
        std::vector<std::vector<std::complex<double>>> temp_hamiltoniano = identidad;
        for (int j = 0; j < N; ++j) {
            if (j == i || j == i + 1) {
                temp_hamiltoniano = kronecker_product(temp_hamiltoniano, pauli_z);
            } else {
                temp_hamiltoniano = kronecker_product(temp_hamiltoniano, identidad);
            }
        }

        #pragma omp critical
        for (int k = 0; k < (1 << N); ++k) {
            for (int l = 0; l < (1 << N); ++l) {
                hamiltoniano[k][l] -= J * temp_hamiltoniano[k][l];
            }
        }
    }

    // Parte relacionada con Pauli-X
    #pragma omp parallel for
    for (int i = 0; i < N; ++i) {
        std::vector<std::vector<std::complex<double>>> temp_hamiltoniano = identidad;
        for (int j = 0; j < N; ++j) {
            if (j == i) {
                temp_hamiltoniano = kronecker_product(temp_hamiltoniano, pauli_x);
            } else {
                temp_hamiltoniano = kronecker_product(temp_hamiltoniano, identidad);
            }
        }

        #pragma omp critical
        for (int k = 0; k < (1 << N); ++k) {
            for (int l = 0; l < (1 << N); ++l) {
                hamiltoniano[k][l] -= G * temp_hamiltoniano[k][l];
            }
        }
    }
}

// Producto de Kronecker con paralelización
std::vector<std::vector<std::complex<double>>> Hamiltoniano::kronecker_product(
    const std::vector<std::vector<std::complex<double>>>& A,
    const std::vector<std::vector<std::complex<double>>>& B) const {
    int rows_A = A.size();
    int cols_A = A[0].size();
    int rows_B = B.size();
    int cols_B = B[0].size();
    std::vector<std::vector<std::complex<double>>> result(
        rows_A * rows_B, std::vector<std::complex<double>>(cols_A * cols_B));

    #pragma omp parallel for collapse(2)
    for (int i = 0; i < rows_A; ++i) {
        for (int j = 0; j < cols_A; ++j) {
            for (int k = 0; k < rows_B; ++k) {
                for (int l = 0; l < cols_B; ++l) {
                    result[i * rows_B + k][j * cols_B + l] = A[i][j] * B[k][l];
                }
            }
        }
    }
    return result;
}
