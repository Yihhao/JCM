from package import *
from package.JCM import JCM_Hamiltonian
from package.plot_energy import plot_energy_b, plot_energy_d, plot_energy_e


def compute(N, walist, wc, g, use_rwa):
    evals_mat = np.zeros((len(walist), N * 2))
    for index, wa in enumerate(walist):
        # evaluate the Hamiltonian
        H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
        # find the energy eigenvalues of the composite system
        evals, ekets = LA.eigh(H)
        evals_mat[index, :] = np.real(evals)
    return evals_mat


def main(N, wc, wa, g, use_rwa, plot_range, savefig):
    wa_range = (plot_range * g) / (2 * pi)
    wc_f = wc / (2 * pi)
    walist = np.linspace(wc_f - wa_range, wc_f + wa_range, 200) * 2 * pi  # atom 1 frequency range
    delta = walist - wc
    print('coupling strength g:%f' % g)

    evals_mat = compute(N, walist, wc, g, use_rwa)

    plot_energy_b(delta, evals_mat, g, wc_f, plot_range, plot_line=5, figsize=(6, 12))
    plt.tight_layout()
    if savefig:
        plt.savefig('./fig/energy_level_b', dpi=720)
    plt.show()

    plot_energy_d(evals_mat, delta, wc, g, plot_range)
    plt.legend(loc=2, prop={'size': 16})
    plt.tight_layout()
    if savefig:
        plt.savefig('./fig/energy_level_d', dpi=720)
    plt.show()

    plot_energy_e(evals_mat, walist, wc, g, delta, plot_range, plot_line=7)
    # plt.legend()
    plt.tight_layout()
    if savefig:
        plt.savefig('./fig/energy_level_e.png')
    plt.show()


if __name__ == '__main__':
    '''
    ref :'PhysRevLett.114.233601'
    plot Fig 1.b
    plot Fig 1.d
    plot Fig 1.e
    '''
    # parameters
    N = 20  # number of cavity fock states
    wc = 1.0 * 2 * pi  # cavity frequency
    wa = 0.0 * 2 * pi  # atom frequency
    g = 0.05 * 2 * pi  # coupling strength
    use_rwa = False  # rwa: rotating wave approximation
    plot_range = 4
    savefig = False
    plot_line = 5
    main(N, wc, wa, g, use_rwa, plot_range, savefig)

    # H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
    # print((wa-wc) / g)
    # evals, ekets = LA.eigh(H)
    # print(np.real(evals[:5]))



