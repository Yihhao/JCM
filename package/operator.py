from numpy import sqrt, zeros, array, arange, tensordot


def density_mtrix(psi):
    psi = psi.reshape(len(psi), 1)
    d = tensordot(psi, psi, axes=([1], [1]))
    return d


def dagger(matrix):
    return matrix.conj().T


def destroy_op(N):
    """annihilation operator for SHO"""
    a = zeros([N, N], dtype=float)
    for i in arange(1, N):
        a[i - 1, i] = sqrt(i)
    return a


def sigmaz_op():
    """Pauli spin 1/2 sigma-z operator."""
    z = array([[1, 0],
               [0, -1]], dtype=float)
    return z


def sigmax_op():
    """Pauli spin 1/2 sigma-x operator."""
    x = array([[0, 1],
               [1, 0]], dtype=float)
    return x


def sigmay_op():
    """Pauli spin 1/2 sigma-y operator."""
    y = array([[0, -1.0j],
               [1.0j, 0]], dtype=complex)
    return y


def sigmap_op():
    """Creation operator for Pauli spins."""
    p = array([[0, 1],
               [0, 0]], dtype=float)
    return p


def sigmam_op():
    """annihilation operator for Pauli spins."""
    m = array([[0, 0],
               [1, 0]], dtype=float)
    return m
