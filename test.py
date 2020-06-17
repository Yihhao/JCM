from package import *
from package.operator import sigmaz_op
from package.JCM import JCM_Hamiltonian, tensor_operator, initial_state
from package.state import fock_state
from package.time_evolution import H_evolution, expect_value


N = 20
# z = 1.414  # 1.7
z = sqrt(25)
# z = 1
wc = 1.0 * 2 * pi
wa = 1.0 * 2 * pi
g = 0.05 * 2 * pi
wav = (0, 1)
use_rwa = False

# psi_text = 'fock'
psi_text = 'coherent'

tlist = np.linspace(0, 50, 5001)

H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
sm, sx, a, I = tensor_operator(N)
psi0, fig_tilte = initial_state(psi_text, N, z, wav)

result = H_evolution(H, psi0, tlist)
op = [dot(sm, dagger(sm)), dot(dagger(sm), sm), dot(dagger(a), a)]
fig, axes = plt.subplots(1, 1, figsize=(12, 8))
ept_g, ept_e, na = expect_value(op, result)
# ept_g = expect_value(sz, result)
# ept_g = expect_value(sz, p_res)
# plt.plot(tlist, np.real(ept_e), label='exit')
plt.plot(tlist, np.real(ept_g))

plt.plot(tlist, np.real(ept_g), 'r', label=r"$|\alpha,e>$")
plt.plot(tlist, np.real(ept_e), 'b', label=r"$|\alpha^{'},g>$")
# plt.plot(tlist, np.real(ept_g), 'r', label=r"$|0,e>$")
# plt.plot(tlist, np.real(ept_e), 'b', label=r"$|1,g>$")
# plt.plot(tlist, np.real(na), 'b', label=r"$Cavity$")
plt.xlim(0, 30)
plt.ylim(0, 1.05)
plt.legend(loc=1, fontsize=20)
plt.xlabel('time', fontsize=20)
plt.ylabel('Probability', fontsize=20)
plt.tight_layout()
# plt.savefig('./fig/rabi_pur_d0.png', dpi=720)
plt.savefig('./fig/rabi_coh_d1.png', dpi=720)
plt.show()
#

fig, axes = plt.subplots(1, 1, figsize=(12, 8))
sz = dot(dagger(sm), sm) - dot(sm, dagger(sm))
ept_g = expect_value(sz, result)
plt.plot(tlist, np.real(ept_g))
plt.xlim(0, 30)
plt.ylim(-1.05, 1.05)
# plt.legend(loc=1)
plt.xlabel('t', fontsize=16)
plt.ylabel(r'$S_z$', fontsize=16)
plt.tight_layout()
# plt.savefig('./fig/rabi_d1.png', dpi=720)
plt.show()