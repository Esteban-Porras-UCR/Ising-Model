#include "Hamiltoniano.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>

// Constructor
Hamiltoniano::Hamiltoniano(int n, double j, double g)
{
    this->N = n;  // Número de spins
    this->J = j;  // Interacción entre spins
    this->G = g;  // Campo magnético
    filas = std::pow(2, N);  // 2^N filas
    columnas = std::pow(2, N);  // 2^N columnas
    matriz_hamiltoniano = std::vector<std::complex<double>>(filas * columnas, {0.0, 0.0});  // Inicializa el Hamiltoniano con ceros

    // Definir las matrices de Pauli y la identidad
    pauli_z = {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), 
               std::complex<double>(0.0, 0.0), std::complex<double>(-1.0, 0.0)}; // Pauli Z
    pauli_x = {std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), 
               std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)}; // Pauli X
    identidad = {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), 
                 std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)}; // Matriz identidad
}

// Método para calcular el producto Kronecker o tensorial de dos matrices
std::vector<std::complex<double>> Hamiltoniano::productoKronecker(const std::vector<std::complex<double>>& A, 
                                                                  const std::vector<std::complex<double>>& B, 
                                                                  int filsA, int colsA, int filsB, int colsB)
{
    std::vector<std::complex<double>> resultado(filsA * filsB * colsA * colsB, {0.0, 0.0});  // Inicializa el resultado con ceros
    for (int i = 0; i < filsA; ++i) {
        for (int j = 0; j < colsA; ++j) {
            for (int k = 0; k < filsB; ++k) {
                for (int l = 0; l < colsB; ++l) {
                    // Realiza el producto Kronecker o tensorial
                    resultado[(i * filsB + k) * (colsA * colsB) + (j * colsB + l)] = A[i * colsA + j] * B[k * colsB + l];
                }
            }
        }
    }
    return resultado;
}

// Método para calcular el producto directo de varias matrices usando el producto Kronecker
std::vector<std::complex<double>> Hamiltoniano::productodirecto(const std::vector<std::vector<std::complex<double>>>& lista_matrices)
{
    std::vector<std::complex<double>> resultado = lista_matrices[0];  // Inicializa el resultado con la primera matriz

    for (int k = 1; k < lista_matrices.size(); ++k) 
    {
        int filasResult = std::sqrt(resultado.size());  // Número de filas de la matriz resultado
        int colsResult = filasResult;  // Número de columnas de la matriz resultado
        int filask = std::sqrt(lista_matrices[k].size());  // Número de filas de la k-ésima matriz
        int colsk = filask;  // Número de columnas de la k-ésima matriz
        // Calcula el producto Kronecker de las matrices
        resultado = productoKronecker(resultado, lista_matrices[k], filasResult, colsResult, filask, colsk);
    }

    return resultado;  // Devuelve el resultado final
}

// Método para construir el Hamiltoniano 
void Hamiltoniano::creacion_Hamiltoniano()
{
    // Parte del Hamiltoniano relacionada con Pauli Z
    for (int i = 0; i < N - 1; i++) 
    {
        std::vector<std::vector<std::complex<double>>> lista_matrices(N, identidad);  // Lista con N matrices identidad
        lista_matrices[i] = pauli_z;  // Sustituye la i-ésima posición por la matriz Pauli Z
        lista_matrices[i + 1] = pauli_z;  // Sustituye la (i+1)-ésima posición por la matriz Pauli Z

        // Calcula el producto directo y actualiza el Hamiltoniano
        std::vector<std::complex<double>> resultado = productodirecto(lista_matrices);
        for (int f = 0; f < filas; ++f) {
            for (int c = 0; c < columnas; ++c) {
                matriz_hamiltoniano[f * columnas + c] -= J * resultado[f * columnas + c];  // Actualiza el Hamiltoniano
            }
        }
    }

    // Parte del Hamiltoniano relacionada con Pauli X
    for (int i = 0; i < N; i++) {
        std::vector<std::vector<std::complex<double>>> lista_matrices(N, identidad);  // Lista con N matrices identidad
        lista_matrices[i] = pauli_x;  // Sustituye la i-ésima posición por la matriz Pauli X
        
        // Calcula el producto directo y actualiza el Hamiltoniano
        std::vector<std::complex<double>> resultado = productodirecto(lista_matrices);
        for (int f = 0; f < filas; ++f) {
            for (int c = 0; c < columnas; ++c) {
                matriz_hamiltoniano[f * columnas + c] -= G * resultado[f * columnas + c];  // Actualiza el Hamiltoniano
            }
        }
    }
}

// Método para mostrar el Hamiltoniano
void Hamiltoniano::mostrar_Hamiltoniano() {
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            std::cout << matriz_hamiltoniano[i * columnas + j] << " ";  // Muestra el Hamiltoniana
        }
        std::cout << std::endl;
    }
}

// Devuelve el Hamiltoniano
std::vector<std::complex<double>> Hamiltoniano::get_matriz_hamiltoniana(){
    return matriz_hamiltoniano;
}

// Constructor de copia
Hamiltoniano::Hamiltoniano(const Hamiltoniano &obj)
{
    // Copia los elementos del Hamiltoniano
    for(int i = 0; i < filas*columnas; ++i)
    {
        matriz_hamiltoniano[i] = obj.matriz_hamiltoniano[i]; 
    }
}