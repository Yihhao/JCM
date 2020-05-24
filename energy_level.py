import matplotlib.pyplot as plt
import numpy as np
from package.JCM import JCM_Hamiltonian
from qutip import *
from numpy import pi, sqrt


def compute(N, walist, wc, g, use_rwa):
    evals_mat = np.zeros((len(walist), N * 2))
    for index, wa in enumerate(walist):
        # evaluate the Hamiltonian
        H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)

        # find the energy eigenvalues of the composite system
        evals, ekets = H.eigenstates()

        evals_mat[index, :] = np.real(evals)

    return evals_mat


def main():
    # initial parameters
    N = 20                  # number of cavity fock states
    wc = 1.0 * 2 * pi       # cavity frequency
    wa = 1.0 * 2 * pi       # atom frequency
    g = 0.05 * 2 * pi       # coupling strength
    use_rwa = True          # rwa: rotating wave approximation
    plot_range = 4
    wa_min = (plot_range * g + wc) / (2*pi)
    walist = np.linspace(-1*wa_min, wa_min, 200) * 2 * pi  # atom 1 frequency range

    print('coupling strength g:%f' % g)

    evals_mat = compute(N, walist, wc, g, use_rwa)

    H0, H1 = JCM_Hamiltonian(N, wc, wa, g, use_rwa, tuple=True)

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
    plt.xlim(-1*plot_range, plot_range)
    plt.yticks([wc_f, 2*wc_f],
               [r'$\omega_c$', r'2$\omega_c$'], fontsize=16)
    plt.tight_layout()
    plt.savefig('./fig/energy_level_1')
    plt.show()
    n = 1
    for i in [(n*2)-1, n*2]:
        c = color[i % 2]
        y = (evals_mat[:, i] - evals_mat[:, 0] - n*wc) / g
        plt.plot((walist - wc) / g, y, c)
    plt.ylim(-4.5, 4.5)
    plt.xlim(-1 * plot_range, plot_range)
    plt.ylabel(r'$\omega / g$')
    plt.savefig('./fig/energy_level_2')
    plt.show()

    H0, H1 = JCM_Hamiltonian(N, wc, wa, g, use_rwa, tuple=True)

    H = H0 + g * H1
    eval = H.eigenenergies()

    plot_energy_levels([H0, g*H1],
                       show_ylabels=True, labels=['|g>', '|g,0>'],
                       N=5, figsize=(8, 12))
    plt.yticks([-0.5*wc, 0.5*wc, 1.5*wc],
               ['|0>', r'$\omega_c   $ |1>', r'2$\omega_c   $ |2>'], fontsize=16)
    plt.text(4.5, eval[1]-g, '|1,\N{MINUS SIGN}>',
             fontdict={'size': 16, 'color': 'k'})
    plt.text(4.5, eval[2]+g, '|1,+>',
             fontdict={'size': 16, 'color': 'k'})
    plt.text(4.5, eval[3]-g, '|2,\N{MINUS SIGN}>',
             fontdict={'size': 16, 'color': 'k'})
    plt.text(4.5, eval[4]+g, '|2,+>',
             fontdict={'size': 16, 'color': 'k'})
    plt.text(4.5, eval[2]-1.5*g, r'$\updownarrow$ 2g',
             fontdict={'size': 20, 'color': 'k'})
    plt.text(4.5, eval[4]-1.7*g, r'$\updownarrow 2\sqrt{2}$g',
             fontdict={'size': 20, 'color': 'k'})

    plt.tight_layout()
    plt.savefig('./fig/energy_level_3', dpi=720)
    plt.show()


if __name__ == '__main__':
    main()
