from qutip import *
from package.state import coherent_state
from package.operator import destory_op
import matplotlib.pyplot as plt
import scipy.linalg as LA
from scipy.linalg import expm
from numpy import sqrt, zeros, array, tensordot
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt


def fock_state(N, n):
    state = zeros([N, 1], dtype=float)
    state[n] = 1
    return state


def plot_fock(rho0, fig=None, ax=None, figsize=(8, 6)):
    if fig is None and ax is None:
        fig, ax = plt.subplots(figsize=(8, 6))
    ax.set_xticks(np.arange(N + 1), minor=False)
    return 0


N = 20
z = 1.2
psi0 = coherent(5, z)

# fig, ax = plot_fock_distribution(psi0)
# ax.set_xticks(np.arange(N + 1), minor=False)
# # plt.savefig('./fig/coherent.png')
# plt.show()

psi = fock_state(5, 0)
print(psi * psi)
psi = psi.reshape(5, 1)
# a = np.arange(60.).reshape(3,4,5)
# b = np.arange(24.).reshape(4,3,2)
# c = np.tensordot(a,b, axes=([1,0],[0,1]))
d = tensordot(psi, psi, axes=([1], [1]))
print(d)
# print(d.shape())
# print(c)
print(ket2dm(psi0))
#
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
#
# u = 1.0
# print(u.conjugate())
# print(coherent(5, 1.0j))
#
# print(coherent_state(5, 1.0j))
# print((destory_op(3)+1.0j).conj().T)
# print((destory_op(3)+1.0j))
