#include "Hamiltoniano.hpp"
#include <cmath>
#include <iostream>
#include <vector>
#include <complex>


// Constructor
Hamiltoniano::Hamiltoniano(int n, double j, double g)
{
    this->N = n;
    this->J = j;
    this->G = g;
    filas = std::pow(2, N);  // 2^N filas
    columnas = std::pow(2, N);  // 2^N columnas
    matriz_hamiltoniano = std::vector<std::complex<double>>(filas * columnas, {0.0, 0.0});

    // Definir matrices de Pauli y la identidad como vectores complejos
    pauli_z = {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), 
               std::complex<double>(0.0, 0.0), std::complex<double>(-1.0, 0.0)}; // Pauli Z
    pauli_x = {std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0), 
               std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0)}; // Pauli X
    identidad = {std::complex<double>(1.0, 0.0), std::complex<double>(0.0, 0.0), 
                 std::complex<double>(0.0, 0.0), std::complex<double>(1.0, 0.0)}; // Matriz identidad
}

// Método para calcular el producto Kronecker
std::vector<std::complex<double>> Hamiltoniano::productoKronecker(const std::vector<std::complex<double>>& A, const std::vector<std::complex<double>>& B, int filsA, int colsA, int filsB, int colsB)
{
    std::vector<std::complex<double>> resultado(filsA * filsB * colsA * colsB, {0.0, 0.0});  // Se asegura de que el tamaño del resultado sea correcto
        for (int i = 0; i < filsA; ++i) {
            for (int j = 0; j < colsA; ++j) {
                for (int k = 0; k < filsB; ++k) {
                    for (int l = 0; l < colsB; ++l) {
                        resultado[(i * filsB + k) * (colsA * colsB) + (j * colsB + l)] = A[i * colsA + j] * B[k * colsB + l];
                    }
                }
            }
        }
    return resultado;
}

// Método para realizar el producto directo de matrices (Producto Kronecker de varias matrices)
std::vector<std::complex<double>> Hamiltoniano::productodirecto(const std::vector<std::vector<std::complex<double>>>& lista_matrices)
{
    std::vector<std::complex<double>> resultado = lista_matrices[0];  // Inicializa el resultado con la primera matriz

    for (int k = 1; k < lista_matrices.size(); ++k) 
    {
        int filasResult = std::sqrt(resultado.size());
        int colsResult = filasResult;
        int filask = std::sqrt(lista_matrices[k].size());
        int colsk = filask;
        resultado = productoKronecker(resultado, lista_matrices[k], filasResult, colsResult, filask, colsk);  // Producto Kronecker // Actualiza el resultado
    }

    return resultado;
}

// Método para construir el Hamiltoniano del modelo Ising
void Hamiltoniano::creacion_Hamiltoniano()
{
    // Parte del Hamiltoniano relacionada con pauli_z
    for (int i = 0; i < N - 1; i++) 
    {
        std::vector<std::vector<std::complex<double>>> lista_matrices(N, identidad);  // Crea una lista con N matrices identidad
        lista_matrices[i] = pauli_z;  // Sustituye la i-ésima posición por la matriz pauli_z
        lista_matrices[i + 1] = pauli_z;  // Sustituye la (i+1)-ésima posición por la matriz pauli_z

        // Calcula el producto directo y actualiza el Hamiltoniano
        std::vector<std::complex<double>> resultado = productodirecto(lista_matrices);
        for (int f = 0; f < filas; ++f) {
            for (int c = 0; c < columnas; ++c) {
                matriz_hamiltoniano[f * columnas + c] -= J * resultado[f * columnas + c];  // Actualiza la matriz Hamiltoniana
            }
        }
    }


    // Parte del Hamiltoniano relacionada con el pauli_x
    for (int i = 0; i < N; i++) {
        std::vector<std::vector<std::complex<double>>> lista_matrices(N, identidad);  // Crea otra lista con N matrices identidad
        lista_matrices[i] = pauli_x;  // Sustituye la i-ésima posición por la matriz pauli_x
        
        // Calcula el producto directo y actualiza el Hamiltoniano
        std::vector<std::complex<double>> resultado = productodirecto(lista_matrices);

        for (int f = 0; f < filas; ++f) {
            for (int c = 0; c < columnas; ++c) {
                matriz_hamiltoniano[f * columnas + c] -= G * resultado[f * columnas + c];  // Actualiza la matriz Hamiltoniana
            }  // Actualiza la matriz Hamiltoniana
        }
    }
}

// Método para mostrar el Hamiltoniano
void Hamiltoniano::mostrar_Hamiltoniano() {
    for (int i = 0; i < filas; ++i) {
        for (int j = 0; j < columnas; ++j) {
            std::cout << matriz_hamiltoniano[i * columnas + j] << " ";
        }
        std::cout << std::endl;
    }
}

std::vector<std::complex<double>> Hamiltoniano::get_matriz_hamiltoniana(){
    return matriz_hamiltoniano;
}


// Constructor de copia
Hamiltoniano::Hamiltoniano(const Hamiltoniano &obj)
{
    for(int i = 0; i < filas*columnas; ++i)
    {
        matriz_hamiltoniano[i] = obj.matriz_hamiltoniano[i]; 
    }
}