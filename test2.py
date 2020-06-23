from package import *
from numpy import arctan, inf, sin, cos, sqrt, exp
from package.JCM import JCM_Hamiltonian, tensor_operator, initial_fock_state
from package.time_evolution import H_evolution, expect
import matplotlib.pyplot as plt


def freq(n, g, delta):
    if delta == 0:
        alpha = arctan(inf)
    else:
        alpha = arctan(2 * g * (n + 1) / delta)
    return alpha


def eigen_energy(n, g, delta):
    omega = sqrt(delta**2+(2*g)**2*(n+1))
    return (wc * (n + 0.5) + 0.5 * omega), (wc * (n + 0.5) - 0.5 * omega)


N = 20  # number of cavity fock states
n = 0  # fock occupy number or amplitude of coherent state
wav = (0, 1)  # Atom initial wavefuction e.g. (0, 1) is |0>
wc = 1.0 * 2 * pi  # cavity frequency
wa = 1.0 * 2 * pi  # atom frequency
g = 0.05 * 2 * pi      # coupling strength that is consistent with chi
delta = wa - wc

use_rwa = True  # rwa: rotating wave approximation

tlist = np.linspace(0, 50, 5001)  # time evolution
Ep, Em = eigen_energy(n, g, delta)
alpha = freq(n, g, delta)

psi_e = []
psi_g = []
for t in tlist:
    A = exp(-1.0j*Ep*t)*(cos(alpha/2)**2) + exp(-1.0j*Em*t)*(sin(alpha/2)**2)
    B = (exp(-1.0j*Ep*t) - exp(-1.0j*Em*t)) * cos(alpha/2) * sin(alpha/2)
    psi_e.append(abs(A)**2)
    psi_g.append(abs(B)**2)


psi0 = initial_fock_state(N, n, wav)
sm, sx, a, I = tensor_operator(N)
H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
res = H_evolution(H, psi0, tlist)
excited = dot(dagger(sm), sm)
ground = dot(sm, dagger(sm))
expt_op = [excited, ground]
excited_list, ground_list = expect(expt_op, res)

plt.plot(tlist[0:-1:50], psi_e[0:-1:50], 'sk:', label='exact solution |0, e>')
# plt.scatter(tlist[0:-1:10], psi_g[0:-1:10], '')
plt.plot(tlist, excited_list, 'r-', label='|0, e>')
plt.plot(tlist, ground_list, 'b-', label='|1, g>')
plt.legend(loc=7, fontsize=14)
plt.xlim(0, 30)
plt.ylim(0, 1.05)
plt.ylabel('Probability', fontsize=14)
plt.xlabel('Time', fontsize=14)
plt.tight_layout()
plt.savefig('./fig/pure_rabi')
plt.show()
