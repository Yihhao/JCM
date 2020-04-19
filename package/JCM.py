from qutip import *


def initial_fock_state(N, n=0, wav=(0, 1), number=1):
    """
    default initial state is |0,1> it mean cavity |0> state and atom |1> state.

    :return psi: a tensor of state |N> * |2>

    :param N: Cavity total state
    :param n: occupy of state cavity
    :param wav: a list of two level system. e.g. (0,1) means |1>
    """
    # TLS: two level system
    TLS = (wav[0] * basis(2, 0) + wav[1] * basis(2, 1)).unit()  # initial state  0*|0> + 1*|1>
    fock = basis(N, n)  # fock state for cavity
    psi = tensor(fock, TLS).unit()
    return psi


def initial_coherent_state(N, n, wav=(1, 0), number=1):
    """
    default initial state is |0,0> it mean cavity coherent state and atom |0> state.

    Parameters
    ----------
    N: Cavity total state
    n: occupy of state cavity
    wav: a list of two level system. e.g. (1,0) means |0>

    Returns
    -------
    psi: a tensor of state |N> * |2>

    """

    TLS = wav[0] * basis(2, 0) + wav[1] * basis(2, 1)
    cavity = coherent(N, n)
    psi = tensor(cavity, TLS).unit()
    return psi


def operator(N, number=1):
    """
    define operator

    :return  :tensor of sm, sz, a and I, the dimension is  |N> * |2>
    """
    j = number + 1  # mean number of 1/2-spin energy level. e.g. 1 paticle of 1/2-spin have 2 energy level
    sm = tensor(qeye(N), destroy(j))
    sz = sm.dag() * sm - sm * sm.dag()  # [S^+, S^-]= S^z
    a = tensor(destroy(N), qeye(j))
    I = tensor(qeye(N), qeye(j))

    return sm, sz, a, I


def JCM_Hamiltonian(N, wc, wa, g, use_rwa=True, number=1, eff=False):
    """Jaynes–Cummings model Hamiltonian"""
    sm, sz, a, I = operator(N, number)
    # define Hamiltonian

    if eff is True:
        delta = wa - wc
        xmi = g ** 2 / delta
        H0 = wc * (a.dag() * a + 0.5 * I) + 0.5 * wa * sz
        H1 = xmi * (a.dag() * a + 0.5 * I) * sz
        return H0, H1

    if use_rwa is True:
        H0 = wc * a.dag() * a + 0.5 * wa * sz
        H1 = g * (sm.dag() * a + sm * a.dag())
    else:
        H0 = wc * a.dag() * a + 0.5 * wa * sz
        H1 = g * (sm.dag() + sm) * (a + a.dag())
    return H0, H1


def time_evolution(H, psi, tlist, N, e_ops=None):
    """
    :param e_ops : expection operators
    """
    sm, sz, a, I = operator(N, number=1)
    if e_ops is None:
        e_ops = [sz, sm.dag() * sm, sm * sm.dag(), sm.dag() * sm + sm * sm.dag()]
    result = sesolve(H, psi, tlist, e_ops)
    return result
