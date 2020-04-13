import matplotlib.pyplot as plt
import numpy as np
# make qutip available in the rest of the notebook
from numpy import pi, sqrt, array
from numpy.linalg import  eigvals
from qutip import *
from JCM import JCM_Hamiltonian, time_evolution, inital_fock_state


def omega_n(delta, g, n=1):
    if delta == 0:
        return 0
    else:
        return sqrt(delta + g**2 * (n+1))

def plot_energies(ng_vec, energies, ymax=(20, 3)):
    """
    Plot energy levels as a function of bias parameter ng_vec.
    """
    fig, axes = plt.subplots(1,1, figsize=(16,6))

    for n in range(len(energies[0,:])):
        axes.plot(ng_vec, energies[:,n])
    axes.set_ylim(-2, ymax[0])
    axes.set_xlabel(r'$n_g$', fontsize=18)
    axes.set_ylabel(r'$E_n$', fontsize=18)

    return fig, axes


if __name__ == '__main__':
    # initial parameters
    wc = 5.4  # cavity frequency
    g = 1e-2  # coupling strength
    N = 2  # number of cavity fock states
    use_rwa = False  # rwa: rotating wave approximation

    # H = JCM_Hamiltonian(N, wc, wa, g)
    deltalist = np.linspace(100, 200*pi, 100) * 1e9
    psi = inital_fock_state(N, n=0)
    sz = np.array([[1,0],
                [0,-1]])
    a = np.array([[0, 1],
                 [0, 0]])
    I = np.eye(2)
    adag = np.array([[0, 0],
                     [1, 0]])
    wa_list = np.linspace(5, 6, 100)
    eigenvalue = []
    evs = []
    for wa in wa_list:
        delta = wa - wc
        H_p = (wc + (g**2 / delta) * sz) * adag * a + 0.5 * (wa + (g**2/ delta)) *sz
        ev = eigvals(H_p)
        eigenvalue.append(ev[1])

    es = []
    for wa in wa_list:
        # wa = wc + delta  # atom frequency
        delta = wa - wc
        H = JCM_Hamiltonian(N, wc, wa, g)
        evals, ekets = H.eigenstates()
        es.append(np.real(evals[1]))

    energies = array([JCM_Hamiltonian(N, wc, wa, g).eigenenergies() for wa in wa_list])
    plot_energies(wa_list, energies)
    plt.show()
    print(energies.shape)
    plt.scatter(wa_list, eigenvalue)
    plt.title('JCM approximation')
    plt.show()
    plt.scatter(wa_list, es)
    plt.title('JCM')
    plt.show()
    energy = []
    energy2 = []
    for wa in wa_list:
        delta = wa - wc
        energy.append(wc*(1+0.5)-0.5*omega_n(delta, g))
        energy2.append(wc * (1 + 0.5) + 0.5 * omega_n(delta, g))
    plt.scatter(wa_list, energy)
    plt.scatter(wa_list, energy2)
    plt.title('exact solution')
    plt.show()