from Jaynes_cummings import Jaynes_cummings
import numpy as np
import matplotlib.pyplot as plt

if __name__=='__main__':
    hbar = 1
    N = 1000
    n = np.array([0, 1, 5, 10, 15, 50, 100])
    t = np.linspace(0, 20, 1000)
    c1 = 0
    c2 = 1
    g = 0.05 * np.pi
    hbar = 1
    wc = 2 * np.pi
    wa = 2 * np.pi
    jaynesC = Jaynes_cummings(N, n, t, c1, c2, g, hbar, wc, wa)
    jaynesC.solver()  # Para QuTiP
    jaynesC.solver_diferenciales()  # Para ecuaciones diferenciales