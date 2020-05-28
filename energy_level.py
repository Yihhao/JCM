import matplotlib.pyplot as plt
import numpy as np
from package.JCM import JCM_Hamiltonian
from package.plot import plot_dispective
from numpy import pi


def compute(N, walist, wc, g, use_rwa):
    evals_mat = np.zeros((len(walist), N * 2))
    for index, wa in enumerate(walist):
        # evaluate the Hamiltonian
        H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()

        evals_mat[index, :] = np.real(evals)

    return evals_mat


def main(N, wc, wa, g, use_rwa, plot_range, savefig):
    wa_min = (plot_range * g + wc) / (2 * pi)
    walist = np.linspace(-1 * wa_min, wa_min, 200) * 2 * pi  # atom 1 frequency range

    print('coupling strength g:%f' % g)

    evals_mat = compute(N, walist, wc, g, use_rwa)

    fig, ax = plt.subplots(figsize=(8, 12))
    color = ['b', 'r', 'k']
    for n in np.arange(5):
        if n == 0:
            c = color[2]
        else:
            c = color[n % 2]
        ax.plot((walist - wc) / g, (evals_mat[:, n] - evals_mat[:, 0]) / (2 * pi), c)
    wc_f = wc / (2 * pi)
    ax.set_xlabel('$\Delta$ /g ', fontsize=16)
    plt.xlim(-1 * plot_range, plot_range)
    plt.yticks([0, wc_f, 2 * wc_f],
               [0, r'$\omega_c$', r'2$\omega_c$'], fontsize=16)
    plt.tight_layout()
    if savefig:
        plt.savefig('./fig/energy_level_1', dpi=720)
    plt.show()
    n = 1
    for i in [(n * 2) - 1, n * 2]:
        c = color[i % 2]
        y = (evals_mat[:, i] - evals_mat[:, 0] - n * wc) / g
        plt.plot((walist - wc) / g, y, c)
    plt.ylim(-7, 7)
    plt.xlim(-1 * plot_range, plot_range)
    plt.ylabel(r'$\omega / g$')
    plt.xlabel(r'$\Delta / g$')
    if savefig:
        plt.savefig('./fig/energy_level_2', dpi=720)
    plt.show()

    H0, Hint = JCM_Hamiltonian(N, wc, wa, g, use_rwa, tuple=True)

    H = H0 + g * Hint
    # H = H0 + Hint
    eval = H0.eigenenergies(eigvals=5)
    print(f'No interatcion:{eval - eval[0]}')
    eval = H.eigenenergies(eigvals=5)
    print(f'interatcion:{eval - eval[0]}')

    plot_line = 5
    plot_dispective([H0, g*Hint], plot_line=plot_line, wa=wa, wc=wc, g=g)
    # plot_dispective([H0, Hint], plot_line=plot_line, wa=wa, wc=wc, g=g)
    plt.tight_layout()
    if savefig:
        plt.savefig('./fig/energy_level3_5.png')
    plt.show()

if __name__ == '__main__':
    # parameters
    N = 20  # number of cavity fock states
    wc = 1.0 * 2 * pi  # cavity frequency
    wa = 2.2 * 2 * pi  # atom frequency
    g = 0.105 * 2 * pi  # coupling strength
    use_rwa = False  # rwa: rotating wave approximation
    plot_range = 7
    savefig = True
    main(N, wc, wa, g, use_rwa, plot_range, savefig)

