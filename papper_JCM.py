import matplotlib.pylab as plt
from package.JCM import JCM_Hamiltonian, initial_coherent_state, operator
from numpy import pi, linspace, array
from qutip import *

N = 20
z = 0
wc = 1.0 * 2 * pi  # cavity freq
wa = 0.5 * 2 * pi
g = 2.0 * 2 * pi  # couple strength
wav = (1, 0)
use_rwa = True
wa_list = array([0, 0.5]) * 2 * pi

psi0 = initial_coherent_state(N, z, wav)
sm, sz, a, I = operator(N)
xc = a.dag() + a
# for wa in wa_list:
tlsit = linspace(0, 10, 1000)
H = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
res = mesolve(H, psi0, tlsit, [], [])
plot_wigner(res.states[0])
plt.show()

xc_list = expect(xc, res.states)
plt.plot(tlsit, xc_list)
plt.show()