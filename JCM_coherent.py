import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
from qutip import *
from package.JCM import JCM_Hamiltonian, operator
from package.state import initial_state


def fig_name():
    # parameter name
    wav_name = ''.join(['_', str(wav[0]), str(wav[1])])
    parameter_name = "".join(['d%s' % int(d), fig_tilte, '%s' % int(z), wav_name])
    detuning = "".join([
        '\n', 'delta= %.2f * 2$\pi$  g =%.2f * 2$\pi$' % (d, couple)
    ])
    # filename = psi_text + parameter_name
    if eff:
        filename = parameter_name + '_eff'
        parameter_name += '_eff'
        text_tilte = "jcm %s = %s" % (fig_tilte, z)
    elif use_rwa:
        filename = parameter_name + '_rwa'
        parameter_name += '_rwa'
        text_tilte = "jcm with rwa %s = %s" % (fig_tilte, z)
    else:
        text_tilte = "jcm without rwa  %s = %s" % (fig_tilte, z)
        filename = parameter_name

    filename = psi_text + filename
    tilte = text_tilte + detuning
    return tilte, filename, parameter_name


def plot():
    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(12, 8))
    fig.suptitle(tilte, fontsize=16)
    axes[0].plot(taulist, nc_list, 'b', label="Cavity")
    axes[0].plot(taulist, na_list, 'r', label="Atom")
    axes[0].set_ylabel("n", fontsize=16)
    axes[0].legend(loc=1)

    axes[1].plot(taulist, na_list, 'r', label="Atom")
    axes[1].set_ylim(0, 1.1)
    axes[1].set_ylabel("n", fontsize=16)
    axes[1].legend(loc=1)

    axes[2].plot(taulist, xc_list, 'b', label="Cavity")
    axes[2].set_ylabel("position", fontsize=16)
    axes[2].legend(loc=1)

    axes[3].plot(taulist, na_list, 'r', label='excited state')
    axes[3].plot(taulist, ground_list, 'b', label="ground state")
    axes[3].set_ylabel('Probability', fontsize=16)
    axes[3].set_xlabel('Time', fontsize=16)
    axes[3].set_ylim(0, 1.1)
    axes[3].legend(loc=1)

    plt.xlim(0, taulist[-1])
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)


if __name__ == '__main__':
    # psi_text = 'fock'
    psi_text = 'coherent'
    store = 'fig'  #store folder
    # initial parameters
    N = 20  # number of cavity fock states
    z = 2  # fock occupy number or amplitude of coherent state
    wav = (1, 1)  # Atom initial wavefuction e.g. (0, 1) is |0>
    wc = 2.0 * 2 * pi  # cavity frequency
    wa = 3.0 * 2 * pi  # atom frequency
    chi = 0.025 * 2 * pi  # parameter in the dispersive hamiltonian >> g**2/delta
    delta = abs(wc - wa)  # detuning
    # g = sqrt(delta * chi)  # coupling strength that is consistent with chi
    g = 0.16 * 2 * pi      # coupling strength that is consistent with chi
    use_rwa = False  # rwa: rotating wave approximation
    eff = True  # eff : use effective Hamiltonian

    taulist = np.linspace(0, 50, 5001)  # time evolution

    # collapse operators
    kappa = 0.005  # cavity dissipation rate
    gamma = 0.05  # atom dissipation rate
    n_th = 0.0  # avg number of thermal bath excitation

    # initial state
    psi0, fig_tilte = initial_state(psi_text, N, z=z, wav=wav)

    print('coupling strength g : %s' % g)
    H = JCM_Hamiltonian(N, wc, wa, g, use_rwa, eff)
    sm, sz, a, I = operator(N)
    c_ops = [sqrt(kappa * (1 + n_th)) * a, sqrt(kappa * n_th) * a.dag()]

    nc = a.dag() * a
    xc = a.dag() + a
    na = sm.dag() * sm
    ground = sm * sm.dag()

    res = mesolve(H, psi0, taulist, [], [])

    nc_list = expect(nc, res.states)
    na_list = expect(na, res.states)
    xc_list = expect(xc, res.states)
    ground_list = expect(ground, res.states)

    d = (delta / (2 * pi))
    couple = (g / (2 * pi))

    path = "".join(['./', store, '/'])

    tilte, filename, parameter_name = fig_name()

    plot()
    # plt.savefig(path+filename, dpi=720)
    plt.show()

    rho_cavity = ptrace(res.states[50], 0)
    plot_wigner(rho_cavity)
    # plt.savefig(path+'winger_'+parameter_name, dpi=720)
    plt.show()

    # plot_wigner(res.states[-1])
    # plt.savefig(path+'psi_winger_'+parameter_name, dpi=720)
    # plt.show()

    rho_cavity_i = ptrace(psi0, 0)
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    plot_fock_distribution(rho_cavity_i, ax=axes[0])
    plot_fock_distribution(rho_cavity, ax=axes[1])
    plt.show()

