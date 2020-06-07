from package import *
# from package.time_evolution import expect_value, H_evolution
from package.operator import sigmaz_op
from package.state import fock_state
from qutip import *


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
    dt_list = np.diff(tlist)
    result.append(rho0)
    for index, dt in enumerate(dt_list):
        # dt = t - dt
        U = expm(-1.0j * H * dt)
        # print('index is %s , error is %s' % (index, LA.norm(U-U_p.full())))
        if two_d:
            temp = dot(dagger(U), rho0)
            temp = dot(temp, U)
            rho0 = temp
            result.appned(temp)
        else:
            temp = dot(U, rho0)
            rho0 = temp
            result.append(temp)
    return result


# hbar = 6.63 * (10 ** (-34)) / (2 * pi)
# me = 9.11 * (10 ** (-31))
# e = 1.6 * (10 ** (-19))
hbar = 1
me = 1
e = 1
B0 = 1
gamma = (e * hbar * B0) / (2 * me)

sz = sigmaz_op()
psi0 = normalization(fock_state(2, 0) + fock_state(2, 1))
H = e / me * B0 * sz
t_list = np.linspace(0, 10, 101)
result = H_evolution(H, psi0, t_list)

sp = np.array([1, 1], dtype=complex) / sqrt(2)
sm = np.array([1, -1], dtype=complex) / sqrt(2)
pp = abs(dot(sp, result)) ** 2
pm = abs(dot(sm, result)) ** 2
plt.plot(t_list, pp)
plt.plot(t_list, pm)
plt.show()
# eigvals, eigvecs = LA.eigh(H)
# psi_t = np.zeros([len(t_list), 2, 2], dtype=complex)

### qutip way
sz = sigmaz()
psi0 = (fock(2, 0) + fock(2, 1)).unit()
H = e / me * B0 * sz
res = sesolve(H, psi0, t_list)

print('time evolution dt = 0: %s ' % LA.norm(result[0] - res.states[0].full()))
print('time evolution dt = 1: %s ' % LA.norm(result[1] - res.states[1].full()))
print('time evolution dt = 2: %s ' % LA.norm(result[2] - res.states[2].full()))
print('time evolution dt = 3: %s ' % LA.norm(result[3] - res.states[3].full()))
print('time evolution dt = 3: %s ' % LA.norm(result[100] - res.states[100].full()))
