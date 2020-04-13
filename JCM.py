from qutip import *


def inital_fock_state(N, n, wav=(0, 1), number=1):
    # TLS: two level system
    TLS = wav[0] * basis(2, 0) + wav[1] * basis(2, 1)  # initial state  0*|0> + 1*|1>
    fock = basis(N, n)  # fock state for cavity
    psi = tensor(fock, TLS).unit()
    return psi


def inital_coherent_state(N, n, wav=(0, 1), number=1):
    TLS = wav[0] * coherent(2, 0) + wav[1] * coherent(2, 1)
    cavity = coherent(N, n)
    psi = tensor(cavity, TLS).unit()
    return psi


def operator(N, number=1):
    """define operator"""
    j = number + 1  # mean number of 1/2-spin energy level. e.g. 1 paticle of 1/2-spin have 2 energy level
    sm = tensor(qeye(N), destroy(j))
    sz = sm.dag() * sm - sm * sm.dag()  # [S^+, S^-]= S^z
    a = tensor(destroy(N), qeye(j))
    return sm, sz, a


def JCM_Hamiltonian(N, wc, wa, g, use_rwa=True, number=1):
    """Jaynes–Cummings model Hamiltonian"""
    sm, sz, a = operator(N, number)
    # define Hamiltonian
    if use_rwa:
        H = wc * a.dag() * a + 0.5 * wa * sz + g * (sm.dag() * a + sm * a.dag())
    else:
        H = wc * a.dag() * a + 0.5 * wa * sz + g * (sm.dag() + sm) * (a + a.dag())
    return H


def time_evolution(H, psi, tlist, N, e_ops=None):
    """e_ops : expection operators"""
    sm, sz, a = operator(N, number=1)
    if e_ops is None:
        e_ops = [sz, sm.dag() * sm, sm * sm.dag(), sm.dag() * sm + sm * sm.dag()]
    result = sesolve(H, psi, tlist, e_ops)
    return result
