from qutip import *
from package.state import coherent_state, fock_state, tensorproduct
from package.operator import destroy_op, dagger, sigmax_op, density_mtrix
import matplotlib.pyplot as plt
import scipy.linalg as LA
from scipy.linalg import expm
from numpy import sqrt, zeros, array, tensordot, arange, real
from package.plot import plot_fock_number
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt


def expct_value(operator, state):
    if len(state) == 1 and len(operator) == 1:
        if operator.shape[0] != operator.shape[0]:
            raise TypeError("input operator dimension must be same [N,N]")
        temp = tensordot(operator, state, axes=[[1], [0]])
        value = tensordot(state, temp, axes=[[0], [0]])
    elif len(operator) == 1:
        value = zeros([len(state)], dtype=complex)
        if operator.shape[0] != operator.shape[0]:
            raise TypeError("input operator dimension must be same [N,N]")
        for idx, s in enumerate(state):
            temp = tensordot(operator, s, axes=[[1], [0]])
            value[idx] = tensordot(s, temp, axes=[[0], [0]])
    else:
        value = zeros([len(operator), len(state)], dtype=complex)
        for iop, op in enumerate(operator):
            if op.shape[0] != op.shape[0]:
                raise TypeError("input operator dimension must be same [N,N]")
            for idx, s in enumerate(state):
                temp = tensordot(op, s, axes=[[1], [0]])
                value[iop, idx] = tensordot(s, temp, axes=[[0], [0]])
    return value


def vn_entropy(rho):
    if len(rho.shape) == 1:
        rho.reshape([len(rho), 1])
    if rho.shape[0] != rho.shape[1]:
        rho = density_mtrix(rho)
        entropy = -1 * (rho * np.log(rho)).trace
    return entropy


# N = 20
# z = 1.0
# psi0 = coherent(5, z)

# fig, ax = plot_fock_distribution(psi0)
# ax.set_xticks(np.arange(N + 1), minor=False)
# # plt.savefig('./fig/coherent.png')
# plt.show()

# plot_wigner(psi0)
# plt.xlabel('X', fontsize=16)
# plt.ylabel('P', fontsize=16)
# plt.savefig('./fig/wigner_coh.png')
# plt.show()
# psi0 = tensor(coherent(N, z), basis(2,1))
# plot_wigner(psi0)
# plt.xlabel('X', fontsize=16)
# plt.ylabel('P', fontsize=16)
# # plt.savefig('./fig/wigner_cohpsi.png')
# plt.show()
# print(destroy_op(5))
N = 5
z = 1.0

a = dagger(destroy_op(N))
psi0 = fock_state(N, 0)
psi1 = fock_state(N, 2)
psi2 = fock_state(N, 1)
state1 = [psi0, psi1, psi2]
# v = tensordot(a, psi0, axes=[[1], [0]])
# print(v)
#
# v = tensordot(psi0, v, axes=[[0], [0]])
# print(v)
v = expct_value(a, state1)
print(v)

a = destroy(N)
# na = a * a.dag()
psi0 = fock(N, 0)
psi1 = fock(N, 2)
psi2 = fock(N, 1)
state2 = [psi0, psi1, psi2]
v = expect(a, state2)
print(v)


