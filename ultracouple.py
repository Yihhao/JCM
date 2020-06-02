# setup the matplotlib graphics library and configure it to show figures inline in the notebook
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
# make qutip available in the rest of the notebook
from qutip import *
from package.JCM import JCM_Hamiltonian
from numpy import pi, sqrt ,real, zeros, array, linspace, where, arange

wc = 1.0 * 2 * pi  # cavity frequency
wa = 1.0 * 2 * pi  # atom frequency

N = 15  # number of cavity fock states
use_rwa = False

# operators
a = tensor(destroy(N), qeye(2))
sm = tensor(qeye(N), destroy(2))

na = sm.dag() * sm  # atom
nc = a.dag() * a  # cavity

# decoupled Hamiltonian
g_vec = linspace(0, 2.0, 101) * 2 * pi  # coupling strength vector

psi_list = []
for index, g in enumerate(g_vec):
    H = JCM_Hamiltonian(N, wc, wa, g, use_rwa=use_rwa)
    # find the groundstate and its energy
    gnd_energy, gnd_state = H.groundstate()

    # store the ground state
    psi_list.append(gnd_state)

na_expt = expect(na, psi_list)  # qubit  occupation probability
nc_expt = expect(nc, psi_list)  # cavity occupation probability

g_idx = where([g_vec == 2 * pi * g for g in [0.0, 0.5, 1.0, 1.5, 2.0]])[1]
psi_sublist = array(psi_list, dtype=object)[g_idx]

xvec = linspace(-5, 5, 200)

fig_grid = (2, len(psi_sublist) * 2)
fig = plt.figure(figsize=(3 * len(psi_sublist), 6))

for idx, psi in enumerate(psi_sublist):
    rho_cavity = ptrace(psi, 0)
    W = wigner(rho_cavity, xvec, xvec)
    ax = plt.subplot2grid(fig_grid, (0, 2 * idx), colspan=2)
    ax.contourf(xvec, xvec, W, 100, norm=mpl.colors.Normalize(-.125, .125), cmap=plt.get_cmap('RdBu'))
    ax.set_title(r"$g = %.1f$" % (g_vec[g_idx][idx] / (2 * pi)), fontsize=16)

# plot the cavity occupation probability in the ground state
ax = plt.subplot2grid(fig_grid, (1, 1), colspan=(fig_grid[1] - 2))
ax.plot(g_vec / (2 * pi), nc_expt, label="Cavity")
ax.plot(g_vec / (2 * pi), na_expt, label="Atom excited state")
ax.legend(loc=0)
ax.set_xlabel('coupling strength(g)')
ax.set_ylabel('Occupation probability')
fig.tight_layout()
plt.savefig('./fig/couple.png', dpi=720)
plt.show()

entropy_cavity = zeros(shape(g_vec))
entropy_atom = zeros(shape(g_vec))

for idx, psi in enumerate(psi_list):
    rho_cavity = ptrace(psi, 0)
    entropy_cavity[idx] = entropy_vn(rho_cavity, 2)

    rho_atom = ptrace(psi, 1)
    entropy_atom[idx] = entropy_vn(rho_atom, 2)

fig, axes = plt.subplots(2, 1, sharex=True ,figsize=(12, 12))
axes[0].plot(g_vec / (2 * pi), nc_expt, 'r', linewidth=2, label="cavity")
axes[0].plot(g_vec / (2 * pi), na_expt, 'b', linewidth=2, label="atom")
axes[0].set_ylabel("Occupation probability", fontsize=16)
# axes[0].set_xlabel("coupling strenght", fontsize=16)
axes[0].legend(loc=0)

axes[1].plot(g_vec/(2*pi), entropy_cavity, 'b', label="cavity", linewidth=2)
axes[1].plot(g_vec/(2*pi), entropy_atom, 'r--', label="atom", linewidth=2)
axes[1].set_ylim(0, 1)
axes[1].set_ylabel("entropy", fontsize=16)
axes[1].set_xlabel("coupling strength(g)", fontsize=16)
axes[1].legend(loc=0)
fig.tight_layout()
plt.savefig('./fig/entropy.png', dpi=720)
plt.show()
