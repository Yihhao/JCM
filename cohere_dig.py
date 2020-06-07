from qutip import *
from package import *
from package.state import coherent_state, fock_state
from package.operator import destroy_op, sigmax_op
from package.plot import plot_fock_number
import matplotlib as mpl
import numpy as np


def vn_entropy(rho):
    if len(rho.shape) == 1:
        rho.reshape([len(rho), 1])
    if rho.shape[0] != rho.shape[1]:
        rho = density_matrix(rho)
        entropy = -1 * (rho * np.log(rho)).trace
    return entropy


def compute(N, walist, wc, g, use_rwa):
    evals_mat = np.zeros((len(walist), N * 2))
    for index, wa in enumerate(walist):
        H0 = wc * nc + 0.5 * dot(wa, sz)
        H1 = dot(sx, a + dagger(a))
        H = H0 + g * H1
        # evaluate the Hamiltonian
        # H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
        # find the energy eigenvalues of the composite system
        # evals, ekets =LA.eigh(H) H.eigenstates()
        evals, ekets = LA.eigh(H)

        evals_mat[index, :] = np.real(evals)

    return evals_mat


def H_evolution(H, rho0, tlist):
    """
    H is independent time
    :param H: Hamilton
    :param rho0: initial state
    :param tlist: a list of time
    :return:
    """

    shape = [len(tlist)]
    two_d = False
    for i in rho0.shape:
        shape.append(i)
    if len(rho0.shape) == 1:
        rho0.reshape([len(rho0), 1])
    if rho0.shape[0] == rho0.shape[1]:
        two_d = True
    result = []
    dt = 0
    # I = eye(H.shape[0])
    for index, t in enumerate(tlist):
        dt = t - dt
        U = expm(-1.0j * H * t)
        if two_d:
            temp = dot(dagger(U), rho0)
            temp = dot(temp, U)
            rho0 = temp
            result.appned(temp)
        else:
            temp = dot(U, rho0)
            result.append(temp)
    return result


nc = 1
sx = 1
sz = 1
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
# v = expct_value(a, state1)
# print(v)

a = destroy(N)
# na = a * a.dag()
psi0 = fock(N, 0)
psi1 = fock(N, 2)
psi2 = fock(N, 1)
state2 = [psi0, psi1, psi2]
v = expect(a, state2)
print(v)


