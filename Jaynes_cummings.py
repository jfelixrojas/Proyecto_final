import numpy as np
import matplotlib.pyplot as plt
from qutip import *
from scipy.integrate import solve_ivp

class Jaynes_cummings:
    """
    Clase que simula el modelo de Jaynes-Cummings usando QuTiP y ecuaciones diferenciales

    Atributos:

    N : int
        Numero máximo de fotones
    n : array
        Numero de fotones para cada simulación
    t : array
        Arreglo de tiempo para la simulación
    c1: float
        Coeficiente del estado inicial para el estado fundamental
    c2: float
        Coeficiente del estado inicial para el estado excitado
    g : float
        Constante de acoplamiento
    hbar : float
        Constante reducida de Planck
    wa : float
        Frecuencia del átomo
    wc : float
        Frecuencia del campo

    Métodos:
    
    solver():
        Soluciona el modelo de Jaynes-Cummings usando QuTiP para cada valor de fotones (n) 
        y grafica la inversión de poblaciones
        
    solver_diferenciales():
        Resuelve el modelo de Jaynes-Cummings usando ecuaciones diferenciales y grafica la inversión de poblaciones
    """

    def __init__(self, N, n, t, c1, c2, g, hbar, wa, wc):
        """ 
        Inicializa la clase Jaynes-Cummings con los parámetros dados.
        
        Parámetros:
        ------------------------------------------------
        N : int
            Numero máximo de fotones
        n : array
            Numero de fotones para cada simulación
        t : array
            Arreglo de tiempo para la simulación
        c1: float
            Coeficiente del estado inicial para el estado fundamental
        c2: float
            Coeficiente del estado inicial para el estado excitado
        g : float
            Constante de acoplamiento
        hbar : float
            Constante reducida de Planck
        wa : float
            Frecuencia del átomo
        wc : float
            Frecuencia del campo
        """
        self.N = N  # Número máximo de fotones
        self.n = n  # Numero de fotones
        self.t = t  # Tiempos
        self.c1 = c1  # Coeficiente del estado base |g>
        self.c2 = c2  # Coeficiente del estado base |e>
        self.g = g  # Constante de acoplamiento
        self.hbar = hbar  # Constante de Planck
        self.wa = wa  # Frecuencia del átomo
        self.wc = wc  # Frecuencia del campo

    def solver(self):
        """ 
        Soluciona el modelo de Jaynes-Cummings usando QuTiP para cada valor de fotones (n) y dibuja la 
        inversión de poblaciones.

        Parámetros:
            self : Jaynes_cummings
            Instancia de la clase con los siguientes atributos
            - N : int
                Número máximo de fotones
            - n : array_like
                Numero de fotones en la simulación
            - t : array_like
                Arreglo de tiempo
            - c1 : complex
                Coeficiente para el estado base |g>.
            - c2 : complex
                Coeficiente para el estado excitado |e>.
            - g : float
                Constante de acople
            - hbar : float
                Constante reducida de Planck
            - wa : float
                Frecuencia del átomo
            - wc : float
                Frecuencia del campo

        Retorno:
        --------
        None
            La función no retorna ningun valor, solo grafica la inversión de poblaciones.
        """
        plt.figure(figsize=(10, 6))  # Configuración de la gráfica

        for i in self.n:
            # Estados iniciales para cada parte del sistema
            psi1 = self.c1 * basis(2, 0) + self.c2 * basis(2, 1)  # Estado inicial del átomo
            psi2 = basis(self.N, i)  # Estado inicial del campo

            # Estado inicial del sistema:
            psi0 = tensor(psi1, psi2).unit()

            # Operadores del problema:
            sigma_plus = tensor(create(2), qeye(self.N))  # Operador de estado excitado
            sigma_minus = tensor(destroy(2), qeye(self.N))  # Operador de estado fundamental
            sigma_z = sigma_plus * sigma_minus - sigma_minus * sigma_plus  # Inversión de poblaciones
            a = tensor(qeye(2), destroy(self.N))  # Operador destrucción 
            a_dagger = tensor(qeye(2), create(self.N))  # Operador creación
            H = ((self.hbar*self.wc * a_dagger * a) + (0.5*self.hbar * self.wa * sigma_z) + self.g * (sigma_plus * a + sigma_minus * a_dagger))
            output = sesolve(H, psi0, self.t, sigma_z, options=None)  # Solución de la ecuación de Schrödinger
            solucion = output.expect[0]

            # Graficar resultados para cada n
            plt.plot(self.g*self.t, solucion, '.-', label='n = {}'.format(i))
        
        plt.title('Inversión de poblaciones')
        plt.grid()
        plt.legend()
        plt.xlabel('g*t')
        plt.ylabel('W(t)')
        plt.savefig('./Imagenes/Simulation_outputs/inversion_poblaciones_qutip.png')
        plt.show()
    
    def solver_diferenciales(self):
        """ 
        Resuelve el modelo de Jaynes-Cummings usando ecuaciones diferenciales y grafica la inversión de poblaciones.

        Parámetros:
            self : Jaynes_cummings
            Instancia de la clase con los siguientes atributos
            - N : int
                Número máximo de fotones
            - n : array_like
                Numero de fotones en la simulación
            - t : array_like
                Arreglo de tiempo
            - c1 : complex
                Coeficiente para el estado base |g>.
            - c2 : complex
                Coeficiente para el estado excitado |e>.
            - g : float
                Constante de acople

        Retorno:
        --------
        None
            La función no retorna ningun valor, solo grafica la inversión de poblaciones.
        """
        plt.figure(figsize=(10, 6))

        def sistema(t, y, n):
            dydt = [-1j * self.g * np.sqrt(n + 1) * y[1],
                    -1j * self.g * np.sqrt(n + 1) * y[0]]
            return dydt

        for i in self.n:
            y0 = [self.c2 + 0j, self.c1 + 0j]
            t_span = (self.t[0], self.t[-1])
            
            # Resolver el sistema utilizando solve_ivp
            solucion = solve_ivp(sistema, t_span, y0, method='RK45', t_eval=self.t, args=(i,))

            # Obtener los valores de tiempo y de la solución
            t = solucion.t
            y1 = solucion.y[0]
            y2 = solucion.y[1]

            # Calcular el módulo cuadrado de y1 y y2
            modulo_cuadrado_y1 = np.abs(y1)**2
            modulo_cuadrado_y2 = np.abs(y2)**2

            # Calcular la resta del módulo cuadrado
            diferencia_modulo_cuadrado = modulo_cuadrado_y1 - modulo_cuadrado_y2
            plt.plot(self.g*t, diferencia_modulo_cuadrado, '.-', label='n = {}'.format(i))
        
        plt.title('Inversión de poblaciones solucionando la ecuacion diferencial')
        plt.xlabel('Tiempo $\lambda t$')
        plt.ylabel(r'$|y_1(t)|^2 - |y_2(t)|^2$')
        plt.legend()
        plt.grid(True)
        plt.savefig('./Imagenes/Simulation_outputs/inversion_poblaciones_diferenciales.png')
        plt.show()