#include "Hamiltoniano.hpp"
#include <sys/time.h>

double seconds() //Función para obtener el tiempo 
{
    struct timeval tmp;
    double sec;
    gettimeofday(&tmp, (struct timezone *)0);
    sec = tmp.tv_sec + ((double)tmp.tv_usec) / 1000000.0;
    return sec;
}

std::vector<std::complex<double>> Ec_Schrodinger(const std::vector<std::complex<double>>& H, 
                                                  const std::vector<std::complex<double>>& onda) 
{
    int N = onda.size();  // Tamaño del estado
    std::vector<std::complex<double>> result(N, {0.0, 0.0});  // Vector resultante

    // Realizamos el producto punto entre H y onda (simulando np.dot)
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            result[i] += -std::complex<double>(0.0, 1.0) * H[i * N + j] * onda[j];
        }
    }
    return result;
}

std::vector<std::complex<double>> rk4(std::vector<std::complex<double>> (*func)(const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& onda), 
                                      const std::vector<std::complex<double>>& H, const std::vector<std::complex<double>>& onda, double h) 
{
    std::vector<std::complex<double>> temp(onda.size(), {0.0, 0.0});

    // Calcular k1
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

    return temp;
}

int main() {
    int N = 10;   // Número de spins
    double J = 1.0; // Interacción entre spins
    double G = 1.0; // Campo magnético

    double time_1 = seconds();
    Hamiltoniano hamilton(N, J, G); // Crea una instancia del Hamiltoniano
    hamilton.creacion_Hamiltoniano(); // Construye el Hamiltoniano
    double time_2 = seconds();
    
    std::vector<std::complex<double>> hamil = hamilton.get_matriz_hamiltoniana();

    int n = 1000;
    double start = 0.0;
    double end = 10.0;

    double h = (end - start) / n;  // Espaciamiento entre puntos


    std::vector<std::complex<double>> onda(std::pow(2, N), {0.0, 0.0});
    onda[0] = std::complex<double>(1.0, 0.0);  // Inicialización de la onda

    // Iteración en el tiempo
    for (int i = 0; i < n; ++i) {
        std::cout << i * h << " " << onda[0].real() << std::endl;  // Imprimir la parte real de la onda
        onda = rk4(Ec_Schrodinger, hamil, onda, h);
    }
    //std::cout << "Tiempo de creación del Hamiltoniano: " << time_2 -time_1 << " segundos" << std::endl;
    return 0;
}