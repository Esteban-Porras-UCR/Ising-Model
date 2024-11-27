#include <vector>
#include <complex>
#include <iostream>

class Hamiltoniano {
public:
    // Constructor para inicializar el Hamiltoniano
    Hamiltoniano(int N, double J, double G);

    // Método para construir el Hamiltoniano
    void crear_Hamiltoniano();

    // Método para mostrar el Hamiltoniano
    void mostrar_Hamiltoniano() const;

    // Método para calcular el producto de Kronecker
    std::vector<std::vector<std::complex<double>>> kronecker_product(
        const std::vector<std::vector<std::complex<double>>>& A,
        const std::vector<std::vector<std::complex<double>>>& B) const;

private:
    int N;  // Número de spins
    double J;  // Interacción entre spins
    double G;  // Campo magnético
    std::vector<std::vector<std::complex<double>>> pauli_x, pauli_z, identidad; // Matrices de Pauli y matriz identidad
    std::vector<std::vector<std::complex<double>>> hamiltoniano;  // Hamiltoniano calculado

    // Método auxiliar para inicializar las matrices de Pauli y la identidad
    void inicializar_matrices_pauli();
};


