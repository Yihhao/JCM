import numpy as np
import scipy.linalg as la
from scipy.linalg import expm
from numpy import pi, sqrt


def unitary_op(t, *args):
    return expm(-1.0j * H * t / hbar)


# hbar = 6.63 * (10 ** (-34)) / (2 * pi)
# me = 9.11 * (10 ** (-31))
# e = 1.6 * (10 ** (-19))
hbar = 1
me = 1
e = 1
B0 = 1
gamma = (e * hbar * B0) / (2 * me)

sz = np.array([[1, 0],
               [0, -1]])
inital_state = np.array([1, 1], dtype=complex) / sqrt(2)
H = e / me * B0 * sz
t_list = np.linspace(0, 10, 200)
eigvals, eigvecs = la.eig(H)
psi_t = np.zeros([len(t_list), 2, 2], dtype=complex)
for idx, t in enumerate(t_list):
    psi_t[idx] = unitary_op(t) * inital_state

expt_val = np.dot(inital_state, psi_t)
    # np.tensordot(inital_state, psi_t, axes=([0], [1, 0]))
