import matplotlib.pyplot as plt
import numpy as np
from numpy import pi, sqrt
from qutip import *
from package.JCM import JCM_Hamiltonian, initial_fock_state, initial_coherent_state, operator

if __name__ == '__main__':
    # text = 'fock'
    text = 'coherent'
    # initial parameters
    N = 20                  # number of cavity fock states
    z = sqrt(4)             # fock occupy number or amplitude of coherent state
    wav = (1, 1)            # Atom initial wavefuction e.g. (0, 1) is |0>
    wc = 2.0 * 2 * pi       # cavity frequency
    wa = 2.0 * 2 * pi       # atom frequency
    chi = 0.025 * 2 * pi    # parameter in the dispersive hamiltonian >> g**2/delta
    delta = abs(wc - wa)    # detuning
    g = sqrt(delta * chi)   # coupling strength that is consistent with chi
    # g = 0.5 * 2 * pi      # coupling strength that is consistent with chi
    use_rwa = False          # rwa: rotating wave approximation
    eff = False              # eff : use effective Hamiltonian

    taulist = np.linspace(0, 50, 5001)  # time evolution

    # collapse operators
    kappa = 0.005  # cavity dissipation rate
    gamma = 0.05   # atom dissipation rate
    n_th  = 0.0    # avg number of thermal bath excitation

    # initial state
    if text == 'coherent':
        psi0 = initial_coherent_state(N, z, wav)
    elif text == 'fock':
        psi0 = initial_fock_state(N, z, wav)
    else:
        psi0 = 0

    print('coupling strength g : %s' % g)
    H = JCM_Hamiltonian(N, wc, wa, g, use_rwa, eff)
    sm, sz, a, I = operator(N)
    c_ops = [sqrt(kappa * (1 + n_th)) * a, sqrt(kappa * n_th) * a.dag()]

    nc = a.dag() * a
    xc = a.dag() + a
    na = sm.dag() * sm

    res = mesolve(H, psi0, taulist, [], [])

    nc_list = expect(nc, res.states)
    na_list = expect(na, res.states)
    xc_list = expect(xc, res.states)
    sz_list = expect(sz, res.states)

    d = (delta / (2 * pi))
    detuning = '  delta=%.2f * 2$\pi$' % d

    if use_rwa:
        text_tilte = "jcm with rwa z = %s" % z
        filename = "".join(['./fig/', text, '_d%sz%s_rwa' % (int(d), int(z))])
    else:
        text_tilte = "jcm without rwa z = %s" % z
        filename = "".join(['./fig/', text, '_d%sz%s' % (int(d), int(z))])

    fig, axes = plt.subplots(3, 1, sharex=True, figsize=(12, 6))
    fig.suptitle(text_tilte + detuning, fontsize=16)
    axes[0].plot(taulist, nc_list, 'b', label="Cavity")
    axes[0].plot(taulist, na_list, 'r', label="Atom")
    axes[0].set_ylabel("n", fontsize=16)
    axes[0].legend(loc=1)

    axes[1].plot(taulist, na_list, 'r', label="Atom")
    axes[1].set_ylabel("n", fontsize=16)
    axes[1].legend(loc=1)

    axes[2].plot(taulist, xc_list, 'b', label="Cavity")
    axes[2].set_ylabel("position", fontsize=16)
    axes[2].set_xlabel('Time', fontsize=16)
    axes[2].legend(loc=1)

    plt.xlim(0, 50)
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    # plt.savefig(filename, dpi=720)
    plt.show()

    # rho_cavity = ptrace(psi.dag() * res.states[-1], 0)
    # plot_wigner(rho_cavity)
    # plt.show()
