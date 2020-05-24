from package.JCM import JCM_Hamiltonian, operator
from qutip import *
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from numpy import pi, sqrt ,real, zeros, array, linspace, where, arange, exp

def log_base(e):
    return e

# parameters
wc = 1.0 * 2 * pi  # cavity frequency
wa = 1.0 * 2 * pi  # atom frequency

N = 15  # number of cavity fock states
use_rwa = False

g_vec = linspace(0, 2.0, 101) * 2 * pi  # coupling strength vector
psi_list = []

H0, Hint = JCM_Hamiltonian(N, wc, wa, use_rwa=use_rwa)
sm, sz, a, I = operator(N)
for index, g in enumerate(g_vec):
    H = H0 + g * Hint
    gnd_energy, gnd_state, = H.groundstate()
    psi_list.append(gnd_state)

na_expt = expect(sm.dag() * sm, psi_list)  # qubit  occupation probability
nc_expt = expect(a.dag() * a, psi_list)  # cavity occupation probability

fig, axes = plt.subplots(1, 1, sharex=True, figsize=(8, 4))

axes.plot(g_vec / (2 * pi), nc_expt, 'r', linewidth=2, label="cavity")
axes.plot(g_vec / (2 * pi), na_expt, 'b', linewidth=2, label="atom")
axes.set_ylabel("Occupation probability", fontsize=16)
axes.set_xlabel("coupling strenght", fontsize=16)
axes.legend(loc=0)

fig.tight_layout()
plt.show()

entropy_cavity = zeros(shape(g_vec))
entropy_atom = zeros(shape(g_vec))

for idx, psi in enumerate(psi_list):
    rho_cavity = ptrace(psi, 0)
    entropy_cavity[idx] = entropy_vn(rho_cavity, 2)

    rho_atom = ptrace(psi, 1)
    entropy_atom[idx] = entropy_vn(rho_atom, 2)
fig, axes = plt.subplots(1, 1, figsize=(12,6))
axes.plot(g_vec/(2*pi), entropy_cavity, 'b', label="cavity", linewidth=2)
axes.plot(g_vec/(2*pi), entropy_atom, 'r--', label="atom", linewidth=2)
axes.set_ylim(0,1)
axes.set_ylabel("entropy", fontsize=16)
axes.set_xlabel("coupling strength", fontsize=16)
axes.legend(loc=0)
plt.show()
