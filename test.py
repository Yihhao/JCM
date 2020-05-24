import matplotlib.pyplot as plt
import numpy as np
from package.JCM import JCM_Hamiltonian
from qutip import *
from numpy import pi, sqrt


def compute(N, walist, wc, g, use_rwa):
    evals_mat = np.zeros((len(walist), N * 2))
    for index, wa in enumerate(walist):
    #     delta = abs(wc - wa)  # detuning
    #     g = sqrt(delta * chi)  # coupling strength that is consistent with chi
        # evaluate the Hamiltonian
        H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()

        evals_mat[index, :] = np.real(evals)

    return evals_mat


def main():
    # initial parameters
    N = 20                  # number of cavity fock states
    wc = 2.0 * 2 * pi       # cavity frequency
    wa = 2.0 * 2 * pi       # atom frequency
    # chi = 0.1 * 2 * pi      # parameter in the dispersive hamiltonian >> g**2/delta
    delta = abs(wc - wa)    # detuning
    # g = sqrt(delta * chi)   # coupling strength that is consistent with chi
    g = 0.16 * 2 * pi
    use_rwa = True          # rwa: rotating wave approximation
    walist = np.linspace(1.0, 3.0, 200) * 2 * pi  # atom 1 frequency range
    print('coupling strength g:%f' % g)

    evals_mat = compute(N, walist, wc, g, use_rwa)
    fig, ax = plt.subplots(figsize=(12, 6))
    for n in [0, 1, 2, 3]:
        ax.plot((walist-wc)/g, (evals_mat[:, n] - evals_mat[:, 0]) / (2 * pi), 'b')

    ax.set_xlabel('Energy splitting of atom')
    ax.set_ylabel('Eigenenergies')
    ax.set_title('Energy spectrum')
    plt.show()

    H0, H1 = JCM_Hamiltonian(N, wc, wa, g, use_rwa, tuple=True)

    H = H0 + H1
    eval = H.eigenenergies()
    print(eval[0: 5])

    plot_energy_levels([H0, H1], N=5, figsize=(8, 4))
    plt.title("Energy level")
    plt.savefig('./fig/energy_level')
    plt.show()


if __name__ == '__main__':
    main()
