#include <iostream>
#include <vector>
#include <complex>
#include <mpi.h>

using namespace std;

// Definición de las matrices de Pauli y la matriz identidad
const vector<vector<complex<double>>> pauli_x = { {0, 1}, {1, 0} }; // Matriz de Pauli X
const vector<vector<complex<double>>> pauli_y = { {0, complex<double>(0, 1)}, {complex<double>(0, -1), 0} }; // Matriz de Pauli Y
const vector<vector<complex<double>>> pauli_z = { {1, 0}, {0, -1} }; // Matriz de Pauli Z
const vector<vector<complex<double>>> identidad = { {1, 0}, {0, 1} }; // Matriz identidad de 2x2

class Hamiltoniano {
public:
    Hamiltoniano(int n, double j, double g) : N(n), J(j), G(g) {
        int size = 1 << N; // Tamaño de la matriz Hamiltoniana (2^N)
        H.resize(size, vector<complex<double>>(size, complex<double>(0.0, 0.0))); // Inicializa la matriz como ceros
    }

    void creacion_Hamiltoniano(int rank, int size) {
        for (int i = rank; i < N - 1; i += size) {
            actualizaHamiltoniano(i);
        }
        for (int i = rank; i < N; i += size) {
            actualizaHamiltonianoPauliX(i);
        }
    }

    void mostrar_Hamiltoniano() const {
        for (const auto& row : H) {
            for (const auto& elem : row) {
                cout << "(" << elem.real() << "," << elem.imag() << ") ";
            }
            cout << endl;
        }
    }

private:
    int N; // Número de spins (o iteraciones)
    double J; // Interacción entre spins
    double G; // Campo magnético
    vector<vector<complex<double>>> H; // Matriz Hamiltoniana

    vector<vector<complex<double>>> productodirecto(const vector<vector<vector<complex<double>>>>& matrices) {
        vector<vector<complex<double>>> resultado = matrices[0]; // Inicializa el resultado con la primera matriz

        for (size_t k = 1; k < matrices.size(); ++k) {
            resultado = kronecker_product(resultado, matrices[k]); // Calcula el producto de Kronecker
        }
        return resultado;
    }

    void actualizaHamiltoniano(int i) {
        vector<vector<vector<complex<double>>>> lista_matrices(N, identidad); // Crea una lista con N matrices identidad

        lista_matrices[i] = pauli_z; // Sustituye la i-ésima posición por la matriz pauli_z
        lista_matrices[i + 1] = pauli_z; // Sustituye la (i+1)-ésima posición por la matriz pauli_z
        
        auto producto = productodirecto(lista_matrices);
        actualizarHamiltonianoConProducto(producto, -J);
    }

    void actualizaHamiltonianoPauliX(int i) {
        vector<vector<vector<complex<double>>>> lista_matrices(N, identidad); // Crea otra lista con N matrices identidad

        lista_matrices[i] = pauli_x; // Sustituye la i-ésima posición por la matriz pauli_x
        
        auto producto = productodirecto(lista_matrices);
        actualizarHamiltonianoConProducto(producto, -G);
    }

    void actualizarHamiltonianoConProducto(const vector<vector<complex<double>>>& producto, double factor) {
        for (int row = 0; row < H.size(); ++row) {
            for (int col = 0; col < H.size(); ++col) {
                H[row][col] += factor * producto[row][col]; // Actualiza el Hamiltoniano
            }
        }
    }

    vector<vector<complex<double>>> kronecker_product(const vector<vector<complex<double>>>& A,
                                                      const vector<vector<complex<double>>>& B) {
        int rowsA = A.size();
        int colsA = A[0].size();
        int rowsB = B.size();
        int colsB = B[0].size();

        vector<vector<complex<double>>> result(rowsA * rowsB, vector<complex<double>>(colsA * colsB));

        for (int i = 0; i < rowsA; ++i) {
            for (int j = 0; j < colsA; ++j) {
                for (int k = 0; k < rowsB; ++k) {
                    for (int l = 0; l < colsB; ++l) {
                        result[i * rowsB + k][j * colsB + l] = A[i][j] * B[k][l]; // Producto de Kronecker
                    }
                }
            }
        }
        return result;
    }
};

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int N = 10;   // Número de espines
    double J = 1.0;   // Interacción entre espines
    double G = 1.0;   // Campo magnético

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    Hamiltoniano hamilton(N, J, G);   // Crea una instancia de la clase Hamiltoniano
    hamilton.creacion_Hamiltoniano(rank, size);   // Construye el Hamiltoniano en paralelo

    if (rank == 0) {   // Solo el proceso raíz imprime el Hamiltoniano
        hamilton.mostrar_Hamiltoniano();   // Imprime el Hamiltoniano
    }

    MPI_Finalize();
    
    return 0;
}
