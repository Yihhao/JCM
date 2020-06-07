from qutip import *
from package import *
from package.state import coherent_state, fock_state
from package.operator import destroy_op, sigmax_op, sigmaz_op, sigmam_op
from package.plot import plot_fock_number
import matplotlib as mpl
from package.plot_energy import plot_energy_b
from package.JCM import JCM_Hamiltonian, operator
from package.time_evolution import expect_value
# from package.time_evolution import expect_value, H_evolution,


def vn_entropy(rho):
    if len(rho.shape) == 1:
        rho.reshape([len(rho), 1])
    if rho.shape[0] != rho.shape[1]:
        rho = density_matrix(rho)
        entropy = -1 * (rho * np.log(rho)).trace
    return entropy


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
            result.append(temp)
    return result


from scipy import integrate
from numpy import pi, inf, exp, linspace, meshgrid
def wigner_function(rho, range=5):
    xvec = linspace(-range, range, 100)
    x,y = np.meshgrid(xvec, xvec)
    if rho.shape[0] != rho.shape[1]:
        rho0 = density_matrix(rho)

    z = -np.hypot(x, y)

    plt.contourf(x, y, z, 100)
    # plt.colorbar()

    # plt.axhline(0, color='white')
    # plt.axvline(0, color='white')

    plt.show()
    # invexp = lambda y: exp(-2*1.0j*p*y)
    #
    # print(integrate.quad(invexp, 0, inf))
    return 0


N = 5
z = 1.7
wc = 2.0 * 2 * pi
wa = 1.0 * 2 * pi
g = 0.15 * 2 * pi

sm = tensorproduct(eye(N), sigmam_op())
sx = tensorproduct(eye(N), sigmax_op())
sz = tensorproduct(eye(N), sigmaz_op())
a = destroy_op(N)
nc = dot(dagger(a), a)
nc_tensor = tensorproduct(nc, eye(2))
a_tensor = tensorproduct(destroy_op(N), eye(2))

nc = dot(dagger(a), a)
tlist = np.linspace(0, 50, 5001)
H0 = wc * nc_tensor + 0.5 * wa * sz
H1 = dot(sx, a_tensor + dagger(a_tensor))
H = H0 + g * H1

# delta = wa -wc
# unitary = np.a
psi0 = tensorproduct(coherent_state(N, z), fock_state(2, 1))
# psi0 = tensorproduct(fock_state(N, z), fock_state(2, 1))

# result = H_evolution(H, psi0, tlist)
# s = dot(dagger(sm), sm)
# # pp = abs(dot(psi0.reshape(N*2), result)) ** 2
# ept_v = expect_value(nc_tensor, result)
# plt.plot(tlist, ept_v)
# plt.show()

sm_p = tensor(qeye(N), sigmam())
sx_p = tensor(qeye(N), sigmax())
sz_p = tensor(qeye(N), sigmaz())
a_p = tensor(destroy(N), qeye(2))
I_p = tensor(qeye(N), qeye(2))

H0_p = wc * a_p.dag() * a_p + 0.5 * wa * sz_p
H1_p = sx_p * (a_p + a_p.dag())
H_p = H0_p + g * H1_p
U_p = (-1.0j * 0.01 * H_p).expm()

# H_p = JCM_Hamiltonian(N, wc, wa, g, False)
# psi0_p = tensor(fock(N, z), fock(2, 1))
psi0_p = tensor(coherent(N, z), fock(2, 1))
sm, sz, a_p, I = operator(N)

# res = mesolve(H_p, psi0_p, tlist, [], [])
res = sesolve(H_p, psi0_p, tlist, [])
ept_p = expect(a_p.dag() * a_p, res.states)

result = H_evolution(H, psi0, tlist)
ept_v = expect_value(nc_tensor, result)

plt.scatter(tlist, ept_v)
plt.scatter(tlist[1], ept_v[1])
plt.show()

fig, ax = plt.subplots(figsize=(12, 6))
plt.scatter(tlist, ept_p)
plt.scatter(tlist[1], ept_p[1])
plt.show()

print('Hamiltonian: %s' % LA.norm(H - H_p.full()))
print('State: %s' % LA.norm(psi0-psi0_p.full()))

print('time evolution dt = 0: %s ' % LA.norm(result[0] - res.states[0].full()))
print('time evolution dt = 1: %s ' % LA.norm(result[1] - res.states[1].full()))
print('time evolution dt = 2: %s ' % LA.norm(result[2] - res.states[2].full()))
print('time evolution dt = 3: %s ' % LA.norm(result[3] - res.states[3].full()))
print('time evolution dt = 4: %s ' % LA.norm(result[500] - res.states[500].full()))
