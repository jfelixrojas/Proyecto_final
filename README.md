## Modelo de Jaynes-Cummings con QuTiP

Este repositorio contiene una implementación del modelo de Jaynes-Cummings utilizando la librería QuTiP (Quantum Toolbox in Python) y ecuaciones diferenciales. El código proporciona una simulación de la dinámica de un sistema compuesto por un átomo de dos niveles y un campo electromagnético en una cavidad resonante.

### Descripción del Código

El código define una clase "Jaynes_cummings" que simula el modelo de Jaynes-Cummings. A continuación se detalla su funcionamiento:

#### Atributos de la Clase

- N (int): Número máximo de fotones.
- n (array): Número de fotones para cada simulación.
- t (array): Arreglo de tiempos para la simulación.
- c1 (float): Coeficiente del estado inicial para el estado fundamental.
- c2 (float): Coeficiente del estado inicial para el estado excitado.
- g (float): Constante de acoplamiento.
- hbar (float): Constante reducida de Planck.
- wa (float): Frecuencia del átomo.
- wc (float): Frecuencia del campo.

#### Métodos

- solver(): Resuelve el modelo de Jaynes-Cummings usando QuTiP para cada valor de fotones (n) y grafica la inversión de poblaciones.
- solver_diferenciales(): Resuelve el modelo de Jaynes-Cummings usando ecuaciones diferenciales y grafica la inversión de poblaciones.

### Instalación y Requisitos

Para ejecutar este código, necesitarás instalar las siguientes librerías de Python:

- NumPy
- Matplotlib
- QuTiP
- SciPy

Puedes instalar estas librerías usando pip:

bash
pip install numpy matplotlib qutip scipy


### Ejecución

1. Clona este repositorio a tu máquina local:

    bash
    git clone https://github.com/jfelixrojas/Proyecto_final/tree/main
    

2. Navega al directorio del repositorio:

    bash
    cd tu-repositorio

3. Ejecuta el código en Python:

    bash
    python3 execution.py


### Resultados

El método 'solver()' grafica la inversión de poblaciones usando QuTiP y guarda la imagen en './Imagenes/Simulation_outputs/inversion_poblaciones_qutip.png'

El método solver_diferenciales() resuelve el modelo utilizando ecuaciones diferenciales y guarda la gráfica de la inversión de poblaciones en './Imagenes/Simulation_outputs/inversion_poblaciones_diferenciales.png'.

### Contribuciones

Las contribuciones al proyecto son bienvenidas. Por favor, abre un "issue" o un "pull request" en GitHub si deseas colaborar.
