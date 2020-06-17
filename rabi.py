from package import *
from package.operator import sigmaz_op
from package.JCM import JCM_Hamiltonian, tensor_operator, initial_state
from package.state import fock_state, coherent_state
from package.time_evolution import H_evolution, expect_value


N = 20
z = 1
wc = 1.0 * 2 * pi
wa = 1.0 * 2 * pi
g = 0.1 * 2 * pi
wav = (0, 1)
use_rwa = False

# psi_text = 'fock'
psi_text = 'coherent'

tlist = np.linspace(0, 50, 5001)

H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
sm, sx, a, I = tensor_operator(N)
psi0, fig_tilte = initial_state(psi_text, N, z, wav)

result = H_evolution(H, psi0, tlist)
ext = tensorproduct(coherent_state(N, z+1), fock_state(2, 0))
gnd = tensorproduct(coherent_state(N, z), fock_state(2, 1))
res_gnd = []
res_ext = []
for r in result:
    temp = tensordot(ext, r, axes=[0, 1])
    temp = tensordot(temp, ext, axes=[1, 0])
    res_ext.append((abs(temp)).reshape([1]))
    temp = tensordot(gnd, r, axes=[0, 1])
    temp = tensordot(temp, gnd, axes=[1, 0])
    res_gnd.append((abs(temp)).reshape([1]))


fig, axes = plt.subplots(1, 1, figsize=(12, 8))
plt.plot(tlist, np.real(res_ext), label=r'$|n,e>$')
plt.plot(tlist, np.real(res_gnd), label=r'$|n+1,g>$')
plt.xlim(0, 20)
plt.ylim(0, 1.05)
plt.legend(loc=1)
plt.xlabel('t', fontsize=16)
plt.ylabel('Probability', fontsize=16)
plt.tight_layout()
plt.savefig('./fig/rabi.png', dpi=720)
plt.show()
