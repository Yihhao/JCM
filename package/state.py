from package.JCM import initial_coherent_state, initial_fock_state
from package.operator import dagger, destory_op
from numpy import zeros
from scipy.linalg import expm


def initial_state(text, N=None, **kwargs):
    if N is None:
        raise ValueError("N must input data")
    kwargs.setdefault('wav', (0, 1))
    kwargs.setdefault('z', 0)

    wav = kwargs['wav']
    z = kwargs['z']
    if text == 'coherent':
        fig_tilte = 'z'
        psi0 = initial_coherent_state(N, z, wav)
    elif text == 'fock':
        fig_tilte = 'n'
        psi0 = initial_fock_state(N, z, wav)
    else:
        psi0 = 0
        fig_tilte = ''
    return psi0, fig_tilte


def fock_state(N, n):
    state = zeros(N, dtype=float)
    state[n] = 1
    return state


def coherent_state(N, z):
    if z == 0:
        return fock_state(N)
    else:
        a = destory_op(N)
        D = expm(z * dagger(a) - z.conjugate() * a)
        return D * fock_state(N, 0)
