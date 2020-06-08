import matplotlib.pyplot as plt
from qutip import *
from numpy import pi, linspace
from package.JCM import JCM_Hamiltonian, operator, initial_coherent_state


N = 20
z = 1.7
wc = 2.0 * 2 * pi
wa = 1.0 * 2 * pi
g = 0.15 * 2 * pi
wav = (1, 1)
use_rwa = False
tlist = linspace(0, 50, 5001)

H_p = JCM_Hamiltonian(N, wc, wa, g, use_rwa)
psi0_p = initial_coherent_state(N, z, wav)
sm, sz, a_p, I = operator(N)

res = sesolve(H_p, psi0_p, tlist, [])
ept_p = expect(a_p.dag() * a_p, res.states)

fig, ax = plt.subplots(figsize=(12, 6))
plt.plot(tlist, ept_p)
plt.show()
