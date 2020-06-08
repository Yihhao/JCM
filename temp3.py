from JCM_2 import JCM_Hamiltonian, operator, initial_coherent_state, initial_state
from package import *
from package.time_evolution import H_evolution, expect_value


N = 20
z = 1.7
wc = 2.0 * 2 * pi
wa = 1.0 * 2 * pi
g = 0.15 * 2 * pi
wav = (1, 1)
use_rwa = False

# psi_text = 'fock'
psi_text = 'coherent'

tlist = np.linspace(0, 50, 5001)

H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
sm, sx, a, I = operator(N)
nc = dot(dagger(a), a)
# psi0 = initial_coherent_state(N, z, wav)
psi0, fig_tilte = initial_state(psi_text, N, z, wav)

result = H_evolution(H, psi0, tlist)
ept_v = expect_value(nc, result)
plt.plot(tlist, ept_v)
plt.show()
