#include "Hamiltoniano.hpp"
#include <iostream>
#include <omp.h>
#include <vector>
#include <complex>

void Hamiltoniano::mostrar_Hamiltoniano() const {  // Implementación del método para mostrar el Hamiltoniano almacenado.
    for (const auto& row : hamiltoniano_matrix) {   // Itera sobre las filas de la matriz.
        for (const auto& elem : row) {   // Itera sobre los elementos de cada fila.
            std::cout << elem << " ";   // Imprime el elemento en la consola.
        }
        std::cout << "\n";  // Nueva línea al final de cada fila.
    }
}

std::vector<std::vector<std::vector<std::complex<double>>>> Hamiltoniano::productodirecto(  // Implementación del método para calcular el producto directo de una lista de matrices.
    const std::vector<std::vector<std::vector<std::complex<double>>>>& matrices) {
    std::vector<std::vector<std::vector<std::complex<double>>>> result;

    #pragma omp parallel for      // Paraleliza la iteración sobre las matrices con OpenMP.
    for (size_t i = 0; i < matrices.size(); ++i) {
        result.push_back(matrices[i]);   // Añade cada matriz al resultado (este es un ejemplo simplificado).
    }
    return result;  // Devuelve el vector resultante.
}

std::vector<std::vector<std::complex<double>>> Hamiltoniano::actualizarHamiltoniano(  // Implementación del método para actualizar el Hamiltoniano.
    int i,
    const std::vector<std::vector<std::complex<double>>>& producto,
    double factor) {

    #pragma omp parallel for      // Paraleliza la iteración sobre las filas de la matriz producto.
    for (size_t j = 0; j < producto.size(); ++j) {
        for (size_t k = 0; k < producto[j].size(); ++k) {    // Itera sobre las columnas de cada fila.
            hamiltoniano_matrix[j][k] += producto[j][k] * factor;   // Actualiza la matriz del Hamiltoniano usando el producto y el factor.
        }
    }
    return hamiltoniano_matrix;  // Devuelve el Hamiltoniano actualizado.
}

std::vector<std::vector<std::complex<double>>> Hamiltoniano::kronecker_product(  // Implementación del método para calcular el producto de Kronecker.
    const std::vector<std::vector<std::complex<double>>>& A,
    const std::vector<std::vector<std::complex<double>>>& B) {
    size_t rowsA = A.size(), colsA = A[0].size();   // Determina las dimensiones de las matrices A y B.
    size_t rowsB = B.size(), colsB = B[0].size();

    std::vector<std::vector<std::complex<double>>> result(rowsA * rowsB, std::vector<std::complex<double>>(colsA * colsB));  // Crea la matriz resultante con dimensiones adecuadas.

    #pragma omp parallel for collapse(2)  // Paraleliza las iteraciones para calcular el producto de Kronecker.
    for (size_t i = 0; i < rowsA; ++i) {
        for (size_t j = 0; j < colsA; ++j) {
            for (size_t k = 0; k < rowsB; ++k) {
                for (size_t l = 0; l < colsB; ++l) {
                   
                    result[i * rowsB + k][j * colsB + l] = A[i][j] * B[k][l];   // Calcula cada elemento del producto de Kronecker.
                }
            }
        }
    }
    return result;  // Devuelve la matriz resultante.
}
