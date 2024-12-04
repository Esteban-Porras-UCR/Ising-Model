#!/usr/bin/env python
import numpy as np
import time
class Clase_Hamiltoniano:
    def __init__(self, N, J, G, matrices_pauli):
        self.N = N  # Número de spins (o iteraciones)
        self.J = J  # Interacción entre spins
        self.G = G  # Campo magnético
        self.pauli_x, self.pauli_y, self.pauli_z, self.identidad = matrices_pauli
        # Inicializa el Hamiltoniano como una matriz de ceros de tamaño 2^N x 2^N
        self.matriz_hamiltoniano = np.zeros((2**N, 2**N))

    # Método para calcular el producto directo de matrices
    def productodirecto(self, lista_matrices):
        resultado = lista_matrices[0]  # Inicializa el resultado con la primera matriz
        for k in lista_matrices[1:]:  # Itera sobre las matrices restantes
            resultado = np.kron(resultado, k)  # Calcula el producto de Kronecker
        return resultado

    # Método para construir el Hamiltoniano del modelo Ising
    def creacion_Hamiltoniano(self):
        # Parte del Hamiltoniano relacionada con pauli_z
        for i in range(self.N - 1):
            lista_matrices = [self.identidad] * self.N  # Crea una lista con N matrices identidad
            lista_matrices[i] = self.pauli_z  # Sustituye la i-ésima posición por la matriz pauli_z
            lista_matrices[i + 1] = self.pauli_z  # Sustituye la (i+1)-ésima posición por la matriz pauli_z
            # Actualiza el Hamiltoniano
            self.matriz_hamiltoniano -= self.J * self.productodirecto(lista_matrices)

        # Parte del Hamiltoniano relacionada con el pauli_x
        for i in range(self.N):
            lista_matrices = [self.identidad] * self.N  # Crea otra lista con N matrices identidad
            lista_matrices[i] = self.pauli_x  # Sustituye la i-ésima posición por la matriz pauli_x
            # Actualiza el Hamiltoniano
            self.matriz_hamiltoniano -= self.G * self.productodirecto(lista_matrices)

    # Método para mostrar el Hamiltoniano
    def mostrar_Hamiltoniano(self):
        print(self.matriz_hamiltoniano)  # Imprime el Hamiltoniano

# Matrices de Pauli y la matriz identidad
paulix = np.array([[0, 1],
                    [1, 0]])  # Matriz de Pauli X
pauliy = np.array([[0, 1j],
                    [-1j, 0]])  # Matriz de Pauli Y
pauliz = np.array([[1, 0],
                    [0, -1]])  # Matriz de Pauli Z
identidad = np.eye(2)  # Matriz identidad de 2x2

Matrices = [paulix, pauliy, pauliz, identidad] # Lista con las matrices necesarias
# Parámetros del modelo Ising
N = np.array([2, 10]) #Vamos a visualizar la dinámica para dos valores de N
Interaccion_de_espines = 1.0  #Interacción entre espines
Campo_magnetico_externo = 1.0  #Campo magnético externo
# Definición de la ecuación de Schrödinger con h = 1
def Ec_Schrodinger(H, onda):
    return -1j * np.dot(H, onda)

# Método de Runge-Kutta de cuarto orden (RK4)
def rk4(func, H, onda, h):
    k1 = h * func(H, onda)
    k2 = h * func(H, onda + 0.5 * k1)
    k3 = h * func(H, onda + 0.5 * k2)
    k4 = h * func(H, onda + k3)

    return onda + (1/6) * (k1 + 2 * k2 + 2 * k3 + k4)

# Configuración de tiempos y condiciones iniciales
times = np.linspace(0.0, 10.0, 1000) # Tiempo
h = times[1] - times[0] # Espaciamiento
resultado = np.zeros(times.size) # Array con el resultado obtenido

t1 = time.time() # Para medir cuanto dura el método
# N es un array con 2 y 10
for n in range(N.size):
    Hamiltoniano = Clase_Hamiltoniano(N[n], Interaccion_de_espines, Campo_magnetico_externo, Matrices)
    # Crea un Hamiltoniano de tamaño (2**N[n],2**N[n])
    Hamiltoniano.creacion_Hamiltoniano()  # Construye el Hamiltoniano

    onda = np.zeros(2**N[n]) # Onda a estudiar
    onda[0] = 1.0 # Estado inicial de la onda

    # Evolución en el tiempo
    with open(f"resultados_N_{N[n]}.dat", 'w') as f:
        for t in range(times.size):
            resultado[t] = onda[0].real

            # Escribir el tiempo y el resultado en el archivo
            f.write(f"{times[t]} {resultado[t]}\n")

            onda = rk4(Ec_Schrodinger, Hamiltoniano.matriz_hamiltoniano, onda, h)

    print(f"El método RK4 para N = {N[n]} duró: {time.time() - t1:.4f}s")
