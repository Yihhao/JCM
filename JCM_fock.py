import matplotlib.pyplot as plt
import numpy as np
# make qutip available in the rest of the notebook
from numpy import pi
from qutip import *
from JCM import JCM_Hamiltonian, time_evolution, inital_fock_state


if __name__ == '__main__':
    # initial parameters
    wc = 2.0 * 2 * pi  # cavity frequency
    delta = 0
    wa = wc + delta  # atom frequency
    g = wa  # coupling strength
    N = 100  # number of cavity fock states
    n = 0    # fock state occupy number of cavity
    use_rwa = False  # rwa: rotating wave approximation

    kappa = 0.005  # cavity dissipation rate
    gamma = 0.05  # atom dissipation rate
    n_th_a = 0.0  # avg number of thermal bath excitation
    t = np.linspace(0, 5, 10001)  # time evolution

    # initial state
    wav = (0, 1)  # initial state  |1>
    psi = inital_fock_state(N, n, wav)

    H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

    output = time_evolution(H, psi, t, N)
    n_a = output.expect[0]
    n_b = output.expect[1]
    n_c = output.expect[2]
    n_d = output.expect[3]

    if use_rwa:
        text_tilte = "jcm with rwa"
    else:
        text_tilte = "jcm without rwa"

    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    axes.plot(g * t, n_a, label="pop inversion")
    axes.legend(loc=0)
    axes.set_xlabel('Time')
    axes.set_ylabel('probability')
    axes.set_title(text_tilte)
    axes.text(0, 0, 'delta=%.2f' % delta)
    axes.axis([0, 20, -1, 1])
    plt.show()

    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    axes.plot(g * t, n_b, label="excited")
    axes.plot(g * t, n_c, label="ground")
    axes.plot(g * t, n_d, label="norm")
    axes.legend(loc=0)
    axes.set_xlabel('Time')
    axes.set_ylabel('probability')
    axes.set_title(text_tilte)
    axes.axis([0, 20, 0, 1.5])
    plt.show()