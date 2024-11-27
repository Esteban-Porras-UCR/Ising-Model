#include "Hamiltoniano.hpp"
#include <iostream>
#include <vector>
#include <complex>

int main() {
    Hamiltoniano hamiltoniano;  // Crea una instancia de la clase Hamiltoniano.

    std::vector<std::vector<std::complex<double>>> matriz1 = {{1.0, 2.0}, {3.0, 4.0}};      // Inicializa dos matrices de ejemplo con valores complejos.
    std::vector<std::vector<std::complex<double>>> matriz2 = {{5.0, 6.0}, {7.0, 8.0}};

    hamiltoniano.mostrar_Hamiltoniano();  // Muestra el Hamiltoniano inicial (que probablemente esté vacío).

    auto resultado = hamiltoniano.kronecker_product(matriz1, matriz2);      // Calcula el producto de Kronecker de las dos matrices de ejemplo.	

    for (const auto& fila : resultado) {  // Imprime el resultado del producto de Kronecker en la consola.
        for (const auto& elem : fila) {
            std::cout << elem << " ";
        }
        std::cout << "\n";  // Nueva línea al final de cada fila.
    }

    return 0;  // Termina el programa correctamente.
}
