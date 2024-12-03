#ifndef HAMILTONIANO_HPP
#define HAMILTONIANO_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <complex>  // Necesario para trabajar con std::complex

class Hamiltoniano 
{
private:
    int N; // Número de spins
    double J; // Interacción entre spins
    double G; // Campo magnético
    int filas;
    int columnas;
    std::vector<std::complex<double>> matriz_hamiltoniano; // Matriz Hamiltoniana
    std::vector<std::complex<double>> pauli_z; // Pauli Z
    std::vector<std::complex<double>> pauli_x; // Pauli X
    std::vector<std::complex<double>> identidad; // Identidad

public:
    Hamiltoniano(int N, double J, double G);  // Constructor
    void creacion_Hamiltoniano();  // Método para crear el Hamiltoniano
    void mostrar_Hamiltoniano();  // Método para mostrar la matriz Hamiltoniana
    std::vector<std::complex<double>> productoKronecker(const std::vector<std::complex<double>>& A, const std::vector<std::complex<double>>& B, int filsA, int colsA, int filsB, int colsB);  // Producto Kronecker con números complejos
    std::vector<std::complex<double>> productodirecto(const std::vector<std::vector<std::complex<double>>>& lista_matrices);  // Producto directo de matrices con números complejos
    std::vector<std::complex<double>> get_matriz_hamiltoniana();  // Obtener la matriz Hamiltoniana
    Hamiltoniano(const Hamiltoniano &obj);
};

#endif // HAMILTONIANO_HPP
