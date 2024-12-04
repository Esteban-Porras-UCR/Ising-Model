#ifndef HAMILTONIANO_HPP
#define HAMILTONIANO_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <complex> 

class Hamiltoniano 
{
private:
    Hamiltoniano(); // Constructor por defecto privado para ocasionar un error
    int N; // Número de spins
    double J; // Interacción entre spins
    double G; // Campo magnético
    int filas; // Número de filas del Hamiltoniano
    int columnas; // Número de columnas del Hamiltoniano
    std::vector<std::complex<double>> matriz_hamiltoniano; // Matriz del Hamiltoniano
    std::vector<std::complex<double>> pauli_z; // Pauli Z
    std::vector<std::complex<double>> pauli_x; // Pauli X
    std::vector<std::complex<double>> identidad; // Identidad

public:
    Hamiltoniano(int N, double J, double G);  // Constructor
    void creacion_Hamiltoniano();  // Método para crear el Hamiltoniano
    void mostrar_Hamiltoniano();  // Método para mostrar la matriz
    std::vector<std::complex<double>> productoKronecker(const std::vector<std::complex<double>>& A, const std::vector<std::complex<double>>& B, int filsA, int colsA, int filsB, int colsB);  // Producto Kronecker 
    std::vector<std::complex<double>> productodirecto(const std::vector<std::vector<std::complex<double>>>& lista_matrices);  // Producto directo de matrices 
    std::vector<std::complex<double>> get_matriz_hamiltoniana();  // Obtener el Hamiltoniana
    Hamiltoniano(const Hamiltoniano &obj); // Constructor de copia
};

#endif // HAMILTONIANO_HPP
