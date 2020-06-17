from package.operator import destroy_op
from package import *


def fock_state(N, n):
    state = zeros([N, 1], dtype=float)
    state[n] = 1
    return state


def coherent_state(N, z):
    if z == 0:
        return fock_state(N, 0)
    else:
        a = destroy_op(N)
        D = expm(z * dagger(a) - z.conjugate() * a)
        alpha = dot(D, fock_state(N, 0))
        return alpha
