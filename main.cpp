#include "Hamiltoniano.hpp"
#include <sys/time.h>

// Función para obtener el tiempo en segundos
double seconds() 
{
    struct timeval tmp;
    double sec;
    gettimeofday(&tmp, (struct timezone *)0);
    sec = tmp.tv_sec + ((double)tmp.tv_usec) / 1000000.0;
    return sec;
}

// Resuelve la ecuación de Schrödinger usando el Hamiltoniana
std::vector<std::complex<double>> Ec_Schrodinger(const std::vector<std::complex<double>>& H, 
                                                  const std::vector<std::complex<double>>& onda) 
{
    int N = onda.size();  // Tamaño de la onda
    std::vector<std::complex<double>> result(N, {0.0, 0.0});  // Vector resultante

    // -i por el producto punto entre el hamiltoniano y la onda 
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i] += -std::complex<double>(0.0, 1.0) * H[i * N + j] * onda[j];
        }
    }
    return result;
}

// Método de Runge-Kutta de 4ta orden para resolver la ecuación de Schrödinger
std::vector<std::complex<double>> rk4(std::vector<std::complex<double>> (*func)(const std::vector<std::complex<double>>& H, 
                                                                                const std::vector<std::complex<double>>& onda), 
                                                                                const std::vector<std::complex<double>>& H, 
                                                                                const std::vector<std::complex<double>>& onda, 
                                                                                double h) 
{
    std::vector<std::complex<double>> temp(onda.size(), {0.0, 0.0});

    // Calcular k1 Se calcula de esta manera para multiplicar cada entrada por h
    std::vector<std::complex<double>> k1 = func(H, onda);
    for (int i = 0; i < onda.size(); ++i) {
        k1[i] *= h;
    }
    for (int i = 0; i < onda.size(); ++i) {
        temp[i] = onda[i] + 0.5 * k1[i];
    }

    // Calcular k2
    std::vector<std::complex<double>> k2 = func(H, temp);
    for (int i = 0; i < onda.size(); ++i) {
        k2[i] *= h;
    }
    for (int i = 0; i < onda.size(); ++i) {
        temp[i] = onda[i] + 0.5 * k2[i];
    }
    
    // Calcular k3
    std::vector<std::complex<double>> k3 = func(H, temp);
    for (int i = 0; i < onda.size(); ++i) {
        k3[i] *= h;
    }
    for (int i = 0; i < onda.size(); ++i) {
        temp[i] = onda[i] + k3[i];
    }

    // Calcular k4
    std::vector<std::complex<double>> k4 = func(H, temp);
    for (int i = 0; i < onda.size(); ++i) {
        k4[i] *= h;
    }
    for (int i = 0; i < onda.size(); ++i) {
        temp[i] = onda[i] + (1.0 / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
    }

    return temp;  // Retorna el nuevo estado de la onda
}

int main() {
    int N = 10;   // Número de spins
    double J = 1.0; // Interacción entre spins
    double G = 1.0; // Campo magnético

    // Medir el tiempo de creación del Hamiltoniano
    double time_1 = seconds();
    Hamiltoniano hamilton(N, J, G); // Crear la instancia del Hamiltoniano
    hamilton.creacion_Hamiltoniano(); // Construir la matriz Hamiltoniana
    double time_2 = seconds();
    
    // Obtener la matriz Hamiltoniana
    std::vector<std::complex<double>> hamil = hamilton.get_matriz_hamiltoniana();

    // Inicialización de parámetros para la simulación
    int n = 1000;
    double start = 0.0;
    double end = 10.0;
    double h = (end - start) / n;  // Espaciado entre puntos

    // Inicialización del estado de la onda
    std::vector<std::complex<double>> onda(std::pow(2, N), {0.0, 0.0});
    onda[0] = std::complex<double>(1.0, 0.0);  // Establecer el valor inicial de la onda

    // Iteración temporal para resolver la ecuación de Schrödinger
    for (int i = 0; i < n; ++i) {
        std::cout << i * h << " " << onda[0].real() << std::endl;  // Mostrar la parte real de la onda
        onda = rk4(Ec_Schrodinger, hamil, onda, h);  // Método RK4
    }
        //Los resultados se redirigen a un archivo .txt y se grafican con gnuplot
    return 0;
}