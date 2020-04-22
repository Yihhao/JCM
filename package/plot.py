import numpy as np
import matplotlib.pyplot as plt

def plot_fock(N, psi):
    rho = []
    n = np.arange(N)
    psi0 = psi.full()
    for i in n:
        prob = abs(psi0[i, 0])**2
        rho.append(prob)
    plt.bar(n, rho)
    plt.ylabel('Occupation probability')
    plt.xlabel('Fock number')
    plt.xlim(0, )
    plt.ylim(0., 1.0)
    plt.show()