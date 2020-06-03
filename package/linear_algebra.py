from package import *
from numpy import trace


def p_trace(rho, size=None, axes=None):
    """
    :param rho: the state is tensor state.
    :param size: the state initial size
    :param axes: int
        to what axis do trace
    :return:

    Examples
    --------
    psi = tensorproduct(fock_state(2, 0), fock_state(5, 0))
    >>>
    partial_trace(psi, (2,5,1,1), axes=0)

    array([[1., 0.],
       [0., 0.]])
    this result is equal to density_matrix(fock_state(2, 0))
    """
    if size[0] * size[1] != rho.shape[0] and size[0] * size[1] != rho.shape[1]:
        raise ValueError("the input size mismatch rho size")
    if size is None:
        raise ValueError("size must input data")
    if axes is None:
        raise ValueError("axes must input data")
    if axes != 0 and axes != 1:
        raise ValueError("axes must 0 or 1")
    if rho.shape[0] != rho.shape[1]:
        rho = density_matrix(rho)

    sub_axes = 0 if axes == 1 else 1
    matrix = zeros([size[axes], size[axes]])
    m = 0
    for i in arange(0, size[axes] * size[sub_axes], size[sub_axes]):
        n = 0
        for j in arange(0, size[axes] * size[sub_axes], size[sub_axes]):
            # trace_matrix = zeros([size[sub_axes], size[sub_axes]])
            trace_matrix = rho[i:i + size[sub_axes], j:j + size[sub_axes]]
            matrix[m, n] = trace(trace_matrix)
            n += 1
        m += 1

    return matrix


def density_matrix(psi):
    psi = psi.reshape(len(psi), 1)
    d = tensordot(psi, psi, axes=([1], [1]))
    return d


def dagger(matrix):
    return matrix.conj().T


def tensorproduct(A, B):
    Na1 = A.shape[0]
    Na2 = A.shape[1]
    A = A.reshape([Na1, 1, Na2, 1])
    Nb1 = B.shape[0]
    Nb2 = B.shape[1]
    B = B.reshape([Nb1, 1, Nb2, 1])
    matrix = tensordot(A, B, axes=[[1, 3], [1, 3]]).transpose([0, 2, 1, 3])
    return matrix.reshape([Na1*Nb1, Na2*Nb2])


def normalization(state):
    if state.shape[1] != 1:
        raise TypeError("input state dimension must be [N,1]")
    sum = 0
    for i in state:
        sum += i**2
    return state / sqrt(sum)
