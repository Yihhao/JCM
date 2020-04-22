import matplotlib.pyplot as plt
import numpy as np
# make qutip available in the rest of the notebook
from numpy import pi
from qutip import *
from package.JCM import JCM_Hamiltonian, time_evolution, initial_fock_state, operator


if __name__ == '__main__':
    # initial parameters
    wc = 2.0 * 2 * pi  # cavity frequency
    delta = 0
    wa = wc + delta  # atom frequency
    g = wa  # coupling strength
    N = 2  # number of cavity fock states
    n = 0    # fock state occupy number of cavity
    use_rwa = False  # rwa: rotating wave approximation

    kappa = 0.005  # cavity dissipation rate
    gamma = 0.05  # atom dissipation rate
    n_th_a = 0.0  # avg number of thermal bath excitation
    t = np.linspace(0, 10, 10001)  # time evolution

    save = False

    # initial state
    wav = (0, 1)  # initial state  |1>
    psi = initial_fock_state(N, n, wav)

    H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

    output = time_evolution(H, psi, t, N)
    n_a = output.expect[0]
    n_b = output.expect[1]
    n_c = output.expect[2]
    n_d = output.expect[3]

    if use_rwa:
        text_tilte = "jcm with rwa n = %s" % n
        filename = "".join(['./fig/fock_rwa_', '%s' % n])
    else:
        text_tilte = "jcm without rwa n = %s" % n
        filename = "".join(['./fig/fock_', '%s' % n])

    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    axes.plot(g * t, n_a, label="pop inversion")
    axes.legend(loc=0)
    axes.set_xlabel('Time')
    axes.set_ylabel('expectation value')
    axes.set_title(text_tilte)
    axes.text(0, 0, 'delta=%.2f' % delta)
    axes.axis([0, 20, -1, 1])
    if save:
        plt.savefig(filename+'_z')
    plt.show()

    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    axes.plot(g * t, n_b, label="excited")
    axes.plot(g * t, n_c, label="ground")
    axes.plot(g * t, n_d, label="norm")
    axes.legend(loc=0)
    axes.set_xlabel('Time')
    axes.set_ylabel('expectation value')
    axes.set_title(text_tilte)
    axes.axis([0, 20, 0, 1.5])
    if save:
        plt.savefig(filename)
    plt.show()

