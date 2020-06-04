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


