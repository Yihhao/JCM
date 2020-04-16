import matplotlib.pyplot as plt
import numpy as np
# make qutip available in the rest of the notebook
from numpy import pi
from qutip import *
from JCM import JCM_Hamiltonian, time_evolution, inital_coherent_state


if __name__ == '__main__':
    # initial parameters
    wc = 2.0 * 2 * pi  # cavity frequency
    delta = 0
    wa = wc + delta  # atom frequency
    g = wa  # coupling strength
    N = 100  # number of cavity fock states
<<<<<<< HEAD
    n = 0  # fock state occupy number of cavity
=======
    n = 5  # fock state occupy number of cavity
>>>>>>> 2ad988864f178ad10fea6f5c3ccf88551f66336e
    use_rwa = True  # rwa: rotating wave approximation

    t = np.linspace(0, 5, 10001)  # time evolution

    # collapse operators
    kappa = 0.005  # cavity dissipation rate
    gamma = 0.05  # atom dissipation rate
    n_th_a = 0.0  # avg number of thermal bath excitation

    # initial state
    psi = inital_coherent_state(N, n)
    H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

    output = time_evolution(H, psi, t, N)
    n_a = output.expect[0]
    n_b = output.expect[1]
    n_c = output.expect[2]
    n_d = output.expect[3]

    if use_rwa:
        text_tilte = "jcm with rwa z = %s" % n
        filename = "".join(['./fig/coherent_rwa_', '%s' % n])
    else:
        text_tilte = "jcm without rwa z = %s" % n
        filename = "".join(['./fig/coherent_', '%s' % n])

    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    axes.plot(g * t, n_a, label="pop inversion")
    axes.legend(loc=0)
    axes.set_xlabel('Time')
    axes.set_ylabel('probability')
    axes.set_title(text_tilte)
    axes.text(0, 0, 'delta=%.2f' % delta)
    axes.axis([0, 40, -1, 1])
    plt.savefig(filename+'_z')
    plt.show()

    fig, axes = plt.subplots(1, 1, figsize=(10, 6))
    axes.plot(g * t, n_b, label="excited")
    axes.plot(g * t, n_c, label="ground")
    axes.plot(g * t, n_d, label="norm")
    axes.legend(loc=0)
    axes.set_xlabel('Time')
    axes.set_ylabel('probability')
    axes.set_title(text_tilte)
    axes.axis([0, 40, 0, 1.5])
    plt.savefig(filename)
    plt.show()