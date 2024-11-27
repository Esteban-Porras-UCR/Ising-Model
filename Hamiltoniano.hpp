#include <vector>
#include <complex>
#include <iostream>

class Hamiltoniano {  // Definición de la clase Hamiltoniano, que encapsula las operaciones relacionadas con matrices y Hamiltonianos.
public:
    void mostrar_Hamiltoniano() const;   // Método para mostrar el Hamiltoniano almacenado.

    std::vector<std::vector<std::vector<std::complex<double>>>> productodirecto(      // Método para calcular el producto directo de una lista de matrices.
        const std::vector<std::vector<std::vector<std::complex<double>>>>& matrices);

    std::vector<std::vector<std::complex<double>>> actualizarHamiltoniano(   // Método para actualizar el Hamiltoniano usando un producto y un factor escalar.
        int i,  // Índice que puede indicar qué parte del Hamiltoniano actualizar.
        const std::vector<std::vector<std::complex<double>>>& producto,  // Matriz con la que se actualiza.
        double factor);  // Factor escalar para la actualización.

    std::vector<std::vector<std::complex<double>>> kronecker_product(      // Método para calcular el producto de Kronecker entre dos matrices.
        const std::vector<std::vector<std::complex<double>>>& A,
        const std::vector<std::vector<std::complex<double>>>& B);

private:
    std::vector<std::vector<std::complex<double>>> hamiltoniano_matrix;      // Matriz interna que representa el Hamiltoniano.
};



