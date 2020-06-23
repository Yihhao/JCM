from package import *
from .operator import sigmax_op, sigmam_op, destroy_op
from .state import coherent_state, fock_state


def initial_state(text, N=None, z=0, wav=(1, 0)):
    if N is None:
        raise ValueError("N must input data")
    if text == 'coherent':
        fig_tilte = 'z'
        psi0 = initial_coherent_state(N, z, wav)
    elif text == 'fock':
        fig_tilte = 'n'
        psi0 = initial_fock_state(N, z, wav)
    else:
        raise TypeError("text must be 'coherent' or 'fock'")
    return psi0, fig_tilte


def initial_fock_state(N, n=0, wav=(0, 1)):
    """
    default initial state is |0,1> it mean cavity |0> state and atom |1> state.

    Parameters
    ----------
    N:
        Cavity total state
    n:
        occupy of state cavity
    wav:
        a list of two level system. e.g. (0,1) means |1>
    Returns
    -------
    psi:
        a tensor of state |N> * |2>
    """
    # TLS: two level system
    TLS = normalization(wav[0] * fock_state(2, 0) + wav[1] * fock_state(2, 1))  # initial state  0*|0> + 1*|1>
    fock_st = fock_state(N, n)  # fock state for cavity
    psi = tensorproduct(fock_st, TLS)
    return psi


def initial_coherent_state(N, z=0., wav=(1, 0)):
    """
    default initial state is |0,0> it mean cavity coherent state and atom |0> state.

    Parameters
    ----------
    N:
        Cavity total state, it will influence a and a.dag dimension.
    z:
        coherent state of exp(za+za*)
    wav:
        a list of two level system. e.g. (1,0) means |0>

    Returns
    -------
    psi:
        a tensor of state |N> * |2>
    """

    TLS = normalization(wav[0] * fock_state(2, 0) + wav[1] * fock_state(2, 1))
    cavity = coherent_state(N, z)
    psi = tensorproduct(cavity, TLS)
    return psi


def tensor_operator(N):
    """
    define operator.
    return tensor of sm, sx, a and I, the dimension is  |N> * |2>

    sm :
        sigmaminus operator
    sx :
        sigmax operator
    a  :
        destroy operator
    I  :
        unit operator
    """
    sm = tensorproduct(eye(N), sigmam_op())
    sx = tensorproduct(eye(N), sigmax_op())
    a = tensorproduct(destroy_op(N), eye(2))
    I = tensorproduct(eye(N), eye(2))

    return sm, sx, a, I


def JCM_Hamiltonian(N, wc, wa, g=None, use_rwa=True, tuple=False):
    """
    Jaynes–Cummings model Hamiltonian

    Parameters
    ----------
    N:
        Cavity total state
    wc:
        frequence of Cavity
    wa:
        frequence of atom
    g:
        couple strength
    use_rwa:
        use rotation wave?
    tuple:
        return a Hamiltonian or [H0, Hint]?
    """

    sm, sx, a, I = tensor_operator(N)
    sz = dot(dagger(sm), sm) - dot(sm, dagger(sm))
    nc = dot(dagger(a), a)
    # define Hamiltonian

    # if eff is True:
    #     delta = wa - wc
    #     xmi = g ** 2 / delta
    #     H0 = wc * (a.dag() * a + 0.5 * I) + 0.5 * wa * sz
    #     H1 = xmi * (a.dag() * a + 0.5 * I) * sz
    if use_rwa is True:
        H0 = wc * nc + 0.5 * wa * sz
        H1 = dot(dagger(sm), a) + dot(sm, dagger(a))
    else:
        H0 = wc * nc + 0.5 * wa * sz
        H1 = dot(sx, a + dagger(a))

    if tuple is True or g is None:
        return H0, H1
    else:
        H = H0 + g * H1
        return H
