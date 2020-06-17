from package.linear_algebra import tensorproduct, density_matrix
from package.state import fock_state, coherent_state
from package.operator import destroy_op, sigmax_op
from package import *
from numpy import dot, zeros, arange
from scipy.linalg import norm





N = 20
n = 1.2
psi0 = coherent_state(N, n)
# psi_p = coherent(N, n)
# print(norm(psi0-psi_p))
a = destroy_op(N)
# sx = sigmax_op()
# d = tensorproduct(a, sx)
# d_p = tensor(a, sx)
# d = tensorproduct(psi0, psi_p)
# d_p = tensor(psi0, psi_p)

# print(d)
# print(d_p)
# print(norm(d-d_p))

# expect value
na = dot(dagger(a), a)
print(expect(na, psi0))
# e = sum(psi0*psi0)
# b = dot(na, psi0)
# temp = sum(psi0 * b)
#
# print(temp)
