from package.JCM import initial_coherent_state, initial_fock_state
from package.operator import dagger, destroy_op
from numpy import zeros, tensordot
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
    state = zeros([N, 1], dtype=float)
    state[n] = 1
    return state


def coherent_state(N, z):
    if z == 0:
        return fock_state(N)
    else:
        a = destroy_op(N)
        D = expm(z * dagger(a) - z.conjugate() * a)
        alpha = tensordot(D, fock_state(N, 0), axes=[[1], [0]])
        return alpha

def tensorproduct(A, B):
    Na1 = A.shape[0]
    Na2 = A.shape[1]
    A = A.reshape([Na1, 1, Na2, 1])
    Nb1 = B.shape[0]
    Nb2 = B.shape[1]
    B = B.reshape([Nb1, 1, Nb2, 1])
    matrix = tensordot(A, B, axes=[[1, 3], [1, 3]]).transpose([0, 2, 1, 3])
    return matrix
