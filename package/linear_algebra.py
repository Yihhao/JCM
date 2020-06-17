from package import *


def density_matrix(psi0):
    dim = len(psi0)
    d = zeros([dim, dim])
    for i, p in enumerate(psi0):
        for j, s in enumerate(psi0):
            d[i, j] = p * s
    return d


def dagger(matrix):
    return matrix.conj().T


def tensor(A, B):
    Na1 = A.shape[0]
    Na2 = A.shape[1]
    A = A.reshape([Na1, 1, Na2, 1])
    Nb1 = B.shape[0]
    Nb2 = B.shape[1]
    B = B.reshape([Nb1, 1, Nb2, 1])
    matrix = tensordot(A, B, axes=[[1, 3], [1, 3]]).transpose([0, 2, 1, 3])
    return matrix.reshape([Na1*Nb1, Na2*Nb2])


def tensorproduct(A, B):
    Na1 = A.shape[0]
    Na2 = A.shape[1]
    Nb1 = B.shape[0]
    Nb2 = B.shape[1]
    matrix = zeros([Na1 * Nb1, Na2 * Nb2])
    l = 0
    for i in arange(0, Na1*Nb1, 2):
        k = 0
        for j in arange(0, Na2*Nb2, 2):
            matrix[i:i + Nb1, j:j + Nb2] = A[l, k] * B
            k += 1
        l += 1
    return matrix.reshape([Na1*Nb1, Na2*Nb2])


def normalization(state):
    if state.shape[1] != 1:
        raise TypeError("input state dimension must be [N,1]")
    sum = 0
    for i in state:
        sum += i**2
    return state / sqrt(sum)
