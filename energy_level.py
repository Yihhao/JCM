import matplotlib.pyplot as plt
import numpy as np
from package.JCM import JCM_Hamiltonian
from qutip import *
from numpy import pi, sqrt


def compute(N, wa_list, wc, g, use_rwa):
    evals_mat = np.zeros((len(wa_list), N * 2))
    for index, wa in enumerate(wa_list):
        # wa = wc - d
        # evaluate the Hamiltonian
        H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()

        evals_mat[index, :] = np.real(evals)

    return evals_mat

if __name__ == '__main__':
    # initial parameters
    N = 20  # number of cavity fock states
    wc = 1.0 * 2 * pi  # cavity frequency
    wa = 1.0 * 2 * pi  # atom frequency
    # chi = 0.1 * 2 * pi  # parameter in the dispersive hamiltonian >> g**2/delta
    # delta = abs(wc - wa)  # detuning
    # g = sqrt(delta * chi)   # coupling strength that is consistent with chi
    g = 0.5 * 2 * pi  # coupling strength that is consistent with chi
    use_rwa = False  # rwa: rotating wave approximation
    wa_list = np.linspace(0.5, 3.0, 200) * 2 * pi  # atom 1 frequency range
    # delta_list = np.linspace(-2.0, 2.0, 200) * 2 * pi  # define |wa-wc| <-> detunning
    print('coupling strength g:%f' % g)

    evals_mat = compute(N, wa_list, wc, g, use_rwa)
    fig, ax = plt.subplots(figsize=(12, 6))
    for n in [0, 1, 2, 3, 4]:
        ax.plot(wa_list/ (2*pi), (evals_mat[:, n] - evals_mat[:, 0]) / (2 * pi), 'b')
        # ax.plot((delta_list/ g), evals_mat[:, n] / (2 * pi), 'b')

    ax.set_xlabel('Energy splitting of atom')
    ax.set_ylabel('Eigenenergies')
    ax.set_title('Energy spectrum')
    plt.show()

    H0, H1 = JCM_Hamiltonian(N, wc, wa, g, use_rwa, tuple=True)

    H = H0 + g * H1
    eval = H.eigenenergies()
    print(eval[0: 5])

    fig, ax = plot_energy_levels([H0, H1],
                                 show_ylabels=True, N=6, figsize=(8, 12))
    plt.title("Energy level")
    # plt.savefig('./fig/energy_level')
    plt.show()
