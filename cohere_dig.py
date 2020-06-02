from qutip import *
from package.state import coherent_state, fock_state, tensorproduct
from package.operator import destroy_op, dagger, sigmax_op
import matplotlib.pyplot as plt
import scipy.linalg as LA
from scipy.linalg import expm
from numpy import sqrt, zeros, array, tensordot, arange
from package.plot import plot_fock_number
import matplotlib as mpl
import numpy as np
import matplotlib.pyplot as plt



# N = 20
# z = 1.0
# psi0 = coherent(5, z)

# fig, ax = plot_fock_distribution(psi0)
# ax.set_xticks(np.arange(N + 1), minor=False)
# # plt.savefig('./fig/coherent.png')
# plt.show()

# print(coherent_state(5, 3))
# print(coherent(5, 3))
# rho0 = coherent_state(5, 1.0)
# plot_fock_number(rho0)
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
# h0 = tensorproduct(destroy_op(N), sigmax_op()).reshape(2*N, 2*N)
# H = tensor(destroy(N), sigmax())
# dH = LA.norm(h0-H.full())
# print(dH)
s = tensorproduct(coherent_state(N, 1.0), fock_state(2, 1))
sp = tensor(coherent(N, 1.0), fock(2, 1))
print(s.reshape(10, 1))
print(sp)
ds = LA.norm(s.reshape(10, 1)-sp.full())
print(ds)


