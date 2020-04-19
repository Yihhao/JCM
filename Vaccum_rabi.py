import matplotlib.pyplot as plt
import numpy as np
from qutip import *
from package.JCM import JCM_Hamiltonian, operator, initial_fock_state


wc = 1.0  * 2 * np.pi  # cavity frequency
wa = 1.0  * 2 * np.pi  # atom frequency
g  = 0.5 #* 2 * np.pi  # coupling strength
N = 15                 # number of cavity fock states
use_rwa = True

tlist = np.linspace(0, 25, 10001)  # time evolution
psi0 = initial_fock_state(N, n=0)
sm, sz, a, I = operator(N)
H0, H1 = JCM_Hamiltonian(N, wc, wa, g)
H = H0 + H1

output = sesolve(H, psi0, tlist, [a.dag() * a, sm.dag() * sm])

fig, ax = plt.subplots(figsize=(8,5))
ax.plot(tlist, output.expect[0], label="Cavity")
ax.plot(tlist, output.expect[1], label="Atom excited state")

ax.legend()
ax.set_xlabel('Time')
ax.set_ylabel('Occupation probability')
ax.set_title('Vacuum Rabi oscillations')
plt.show()

plt.plot(tlist, (np.cos(tlist*g))**2, 'r--')
plt.show()
