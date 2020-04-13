import numpy as np
import matplotlib.pyplot as plt
from numpy import sqrt, array
from JCM import JCM_Hamiltonian, time_evolution, inital_fock_state
from numpy.linalg import  eigvals

def omega_n(delta, g, n=0):
    if delta == 0:
        return 0
    else:
        return sqrt(delta**2 + g**2 * (n+1))


if __name__ == '__main__':
    # initial parameters
    wc = 5.4  # cavity frequency
    g = 1e-2  # coupling strength
    N = 2  # number of cavity fock states
    use_rwa = False  # rwa: rotating wave approximation

    energy= []
    wa_list = np.linspace(4.5, 6.5, 100)
    energies = array([JCM_Hamiltonian(N, wc, wa, g).eigenenergies() for wa in wa_list])

    fig, axes = plt.subplots(1,1, figsize=(8,6))

    for n in [1, 2]:
        plt.plot(wa_list, energies[:, n])
    axes.set_xlim(wa_list[0], wa_list[-1])
    axes.set_xlabel(r'$w_a$', fontsize=18)
    axes.set_ylabel(r'$E_n$', fontsize=18)
    plt.show()
    """another way plot anti-crossing figure"""
    eigenvalue = []
    sz = np.array([[1, 0],
                   [0, -1]])
    a = np.array([[0, 1],
                  [0, 0]])
    I = np.eye(2)
    adag = np.array([[0, 0],
                     [1, 0]])
    for wa in wa_list:
        delta = wa - wc
        H_p = (wc + (g**2 / delta) * sz) * adag * a + 0.5 * (wa + (g**2/ delta)) *sz
        ev = eigvals(H_p)
        eigenvalue.append(ev[1])
    plt.scatter(wa_list, eigenvalue)
    plt.title('JCM approximation')
    plt.show()