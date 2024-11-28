#include "Hamiltoniano.hpp"
#include <iostream>
#include <vector>
#include <complex>

int main() {
    // Crear un objeto Hamiltoniano con N = 2, J = 1.0, G = 0.5
    Hamiltoniano hamiltoniano(10, 1.0, 1.0);

    // Construir el Hamiltoniano del modelo
    hamiltoniano.crear_Hamiltoniano();

    // Mostrar el Hamiltoniano
    std::cout << "Hamiltoniano generado:\n";
    hamiltoniano.mostrar_Hamiltoniano();

    return 0;  // Finalizar el programa correctamente
}
