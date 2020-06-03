from qutip import *
from package import *
from package.state import coherent_state, fock_state
from package.operator import destroy_op, sigmax_op
from package.plot import plot_fock_number
import matplotlib as mpl
from package.time_evolution import H_evolution


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
        rho = density_matrix(rho)
        entropy = -1 * (rho * np.log(rho)).trace
    return entropy


from numpy import trace





N = 4
z = 1.0
tlist = np.linspace(0, 5, 10)

a1 = destroy_op(N)
x = sigmax_op()
# psi0 = fock_state(N, 0)
psi0 = tensorproduct(fock_state(N, 2), fock_state(2, 1))
# H1 = np.dot(dagger(a1), a1)
# H1 = x
H1 = tensorproduct(a1, x)
res = H_evolution(H1, psi0, tlist)
# print('A : %s' % res[5])

rho = density_matrix(psi0)
from numpy import trace

axes = 0
size = (3, 2, 1, 1)
# sub_axes = 0 if axes == 1 else 1
# # size[axes] = 5
# # size[sub_axes] = 2
# matrix = zeros([size[axes], size[axes]])
# m = 0
# for i in arange(0, size[axes] * size[sub_axes], size[sub_axes]):
#     n = 0
#     for j in arange(0, size[axes] * size[sub_axes], size[sub_axes]):
#         # trace_matrix = zeros([size[sub_axes], size[sub_axes]])
#         trace_matrix = rho[i:i + size[sub_axes], j:j + size[sub_axes]]
#         matrix[m, n] = trace(trace_matrix)
#         n += 1
#     m += 1
# if size[0] * size[1] != rho.shape[0] or size[0] * size[1] != rho.shape[1]:
#     print(size[0] * size[1])
#     print(rho.shape[0])
#     print(rho.shape[1])
    # raise ValueError("the input size mismatch rho size")
print(p_trace(psi0, (4, 2, 1, 1), 0))

a = destroy(N)
x = sigmax()
# psi0 = fock(N, 0)
psi0 = tensor(fock(N, 2), fock(2, 1))
print(ptrace(psi0, 0))
# H2 = a.dag() * a
# H2 = x
# H2 = tensor(a, x)
# # U = (-1.0j * H2 * 2).expm()
# # rho = U.dag() * ket2dm(psi0) * U
# result = mesolve(H2, psi0, tlist, [], [])
# print(result.states[5])
#
# print(LA.norm(H1 - H2.full()))


# psi0 = fock_state(N, 0)
# psi1 = fock_state(N, 2)
# psi2 = fock_state(N, 1)
# state1 = [psi0, psi1, psi2]
# v = tensordot(a, psi0, axes=[[1], [0]])
# print(v)
#
# v = tensordot(psi0, v, axes=[[0], [0]])
# print(v)
# v = expct_value(a, state1)
# print(v)
#
# a = destroy(N)
# # na = a * a.dag()
# psi0 = fock(N, 0)
# psi1 = fock(N, 2)
# psi2 = fock(N, 1)
# state2 = [psi0, psi1, psi2]
# v = expect(a, state2)
# print(v)
